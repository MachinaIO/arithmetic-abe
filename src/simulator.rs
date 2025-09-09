use std::{result, sync::Arc};

use bigdecimal::{BigDecimal, FromPrimitive, One};
pub use mxx::simulator::lattice_estimator::run_lattice_estimator_cli;
use mxx::{
    circuit::PolyCircuit,
    poly::dcrt::poly::DCRTPoly,
    simulator::{
        SimulatorContext,
        lattice_estimator::{Distribution, EstimatorCliError},
        poly_matrix_norm::PolyMatrixNorm,
        wire_norm::NormPltLweEvaluator,
    },
};
use num_bigint::BigUint;
use thiserror::Error;
use tracing::info;

#[derive(Debug, Error)]
pub enum SimulatorError {
    /// Error bubbling up from the lattice estimator CLI.
    #[error(transparent)]
    Estimator(#[from] EstimatorCliError),
    /// No knapsack length in [1, max_knapsack_len] reached the target security.
    #[error(
        "secure ring dimension len not found for secpar={target_secpar}, crt_bits={crt_bits}, crt_depth={crt_depth}, base_bits={base_bits}"
    )]
    RingDimNotFound {
        target_secpar: u32,
        crt_bits: u32,
        crt_depth: u32,
        base_bits: u32,
        log_dim_range: (u32, u32),
    },
    /// No knapsack length in [1, max_knapsack_len] reached the target security.
    #[error(
        "secure knapsack len not found for secpar={target_secpar}, ring_dim={ring_dim}, max_knapsack_len={max_knapsack_len}, q={q}"
    )]
    KnapsackNotFound {
        target_secpar: u32,
        ring_dim: BigUint,
        max_knapsack_len: u32,
        q: BigUint,
    },
    /// Could not find a log_alpha that reaches the target security.
    #[error(
        "good log_alpha not found for target_secpar={target_secpar}, ring_dim={ring_dim}, log_q={log_q}, m={m}"
    )]
    LogAlphaNotFound {
        target_secpar: u32,
        ring_dim: BigUint,
        log_q: u32,
        m: BigUint,
    },
    #[error("correctness does not hold: error={e}, q_over_4={q_over_4}")]
    NotCorrect { e: BigDecimal, q_over_4: BigDecimal },
}

// #[derive(Debug, Clone, Copy)]
// enum ParamChange {
//     Up,
//     Down,
// }

// Output (crt_depth, base_bits, log_dim, e_b_log_alpha, knapsack_len) or None
pub fn bruteforce_params(
    target_secpar: u32,
    crt_bits: u32,
    crt_depth_range: (u32, u32),
    base_bits_range: (u32, u32),
    log_dim_range: (u32, u32),
    circuit: PolyCircuit<DCRTPoly>,
    input_size: usize,
) -> Option<(u32, u32, u32, i64, u32)> {
    // (cost, crt_depth, base_bits, log_dim, e_b_log_alpha, knapsack_len)
    let mut outputs = Vec::<(u32, u32, u32, u32, i64, u32)>::new();
    for base_bits in base_bits_range.0..=base_bits_range.1 {
        let mut lo = crt_depth_range.0;
        let mut hi = crt_depth_range.1;
        while lo <= hi {
            let crt_depth = lo + ((hi - lo) / 2);
            let (log_dim, e_b_log_alpha, knapsack_len) = match find_min_ring_dim(
                target_secpar,
                crt_bits,
                crt_depth,
                base_bits,
                log_dim_range,
            ) {
                Ok(result) => result,
                Err(e) => {
                    info!(
                        "Security error with target_secpar = {target_secpar}, crt_bits = {crt_bits}, base_bits = {base_bits}, crt_depth = {crt_depth}, input_size = {input_size}: {e}"
                    );
                    // try smaller crt_depth
                    hi = crt_depth - 1;
                    continue;
                }
            };
            match check_correctness(
                target_secpar,
                log_dim,
                crt_bits,
                crt_depth,
                base_bits,
                knapsack_len,
                e_b_log_alpha,
                &circuit,
                input_size,
            ) {
                Ok(cost) => {
                    info!(
                        "Found with target_secpar = {target_secpar}, crt_bits = {crt_bits}, base_bits = {base_bits}, crt_depth = {crt_depth}, input_size = {input_size}, cost = {cost}"
                    );
                    outputs.push((
                        cost,
                        crt_depth,
                        base_bits,
                        log_dim,
                        e_b_log_alpha,
                        knapsack_len,
                    ))
                }
                Err(e) => {
                    info!(
                        "Correctness error with target_secpar = {target_secpar}, crt_bits = {crt_bits}, base_bits = {base_bits}, crt_depth = {crt_depth}, input_size = {input_size}: {e}"
                    );
                    // try larger crt_depth
                    lo = crt_depth + 1;
                }
            }
        }
    }
    outputs
        .into_iter()
        .min_by(|x, y| x.0.cmp(&y.0))
        .map(|outs| (outs.1, outs.2, outs.3, outs.4, outs.5))
}

fn find_min_ring_dim(
    target_secpar: u32,
    crt_bits: u32,
    crt_depth: u32,
    base_bits: u32,
    log_dim_range: (u32, u32),
) -> Result<(u32, i64, u32), SimulatorError> {
    for log_dim in log_dim_range.0..=log_dim_range.1 {
        let ring_dim = BigUint::from(2u32).pow(log_dim);
        match check_security(target_secpar, &ring_dim, crt_bits, crt_depth, base_bits) {
            Ok((log_alpha, knapsack_len)) => {
                return Ok((log_dim, log_alpha, knapsack_len));
            }
            Err(_) => {
                continue;
            }
        }
    }
    Err(SimulatorError::RingDimNotFound {
        target_secpar,
        crt_bits,
        crt_depth,
        base_bits,
        log_dim_range,
    })
}

fn check_security(
    target_secpar: u32,
    ring_dim: &BigUint,
    crt_bits: u32,
    crt_depth: u32,
    base_bits: u32,
) -> Result<(i64, u32), SimulatorError> {
    let log_q = (crt_bits * crt_depth) as u32;
    let q = BigUint::from(2u32).pow(log_q as u32);
    let m_g = crt_bits.div_ceil(base_bits) * crt_depth;
    let m_b = m_g + 2;
    // The column size of the matrix B (sampled with a trapdoor) is m_b; however, one column is an identity polynomial, so we need to ignore one column.
    // Additionally, one more uniformly random matrix is used for encrypting a message in ABE; thus the total column size for ring-LWE is m_b - 1 + 1 = m_b.
    let log_alpha =
        find_log_alpha_for_ring_lwe(target_secpar, ring_dim, log_q, &BigUint::from(m_b))?;
    let knapsack_len = find_knapsack_len(target_secpar, ring_dim, &q, m_b - 1)?;
    Ok((log_alpha, knapsack_len))
}

/// Returns the smallest `knapsack_len` in [1, max_knapsack_len] whose estimated
/// security is at least `target_secpar`, or an error if estimation fails or none found.
/// - `target_secpar`: required minimum security parameter.
/// - `ring_dim`: base ring dimension.
/// - `q`: modulus (as BigUint).
/// - `max_knapsack_len`: upper bound to search (inclusive).
fn find_knapsack_len(
    target_secpar: u32,
    ring_dim: &BigUint,
    q: &BigUint,
    max_knapsack_len: u32,
) -> Result<u32, SimulatorError> {
    for knapsack_len in 1..=max_knapsack_len {
        // Effective LWE dimension n = ring_dim * knapsack_len - ring_dim
        let n = ring_dim * BigUint::from(knapsack_len) - ring_dim;
        // s_dist = Ternary, e_dist = Ternary, m = n, exact = false (rough)
        let secpar = run_lattice_estimator_cli(
            &n,
            q,
            &Distribution::Ternary,
            &Distribution::Ternary,
            Some(&n),
            false,
        )?;

        if secpar as u32 >= target_secpar {
            return Ok(knapsack_len);
        }
    }

    Err(SimulatorError::KnapsackNotFound {
        target_secpar,
        ring_dim: ring_dim.clone(),
        max_knapsack_len,
        q: q.clone(),
    })
}

/// Binary-search for the smallest integer `log_alpha` in [-log_q, -1] such that
/// the estimated security for ring-LWE with parameters (ring_dim, q=2^log_q,
/// s_dist=Ternary, e_dist=DiscreteGaussianAlpha(alpha=2^{-log_alpha}), m) is at
/// least `target_secpar`.
///
/// Returns the found `log_alpha` on success, or an error if estimation fails or
/// no such `log_alpha` exists in the search range.
fn find_log_alpha_for_ring_lwe(
    target_secpar: u32,
    ring_dim: &BigUint,
    log_q: u32,
    m: &BigUint,
) -> Result<i64, SimulatorError> {
    // q = 2^{log_q}
    let q = BigUint::from(1u8) << (log_q as usize);

    // Search bounds (inclusive) over integer log_alpha.
    let mut lo: i64 = -(log_q as i64);
    let mut hi: i64 = -1;
    let mut found: Option<i64> = None;

    while lo <= hi {
        let mid = lo + ((hi - lo) / 2);

        // alpha = sigma/q = 2^{log_alpha}
        let alpha = 2f64.powi(mid as i32); // safe for practical parameter sizes

        let e_dist = Distribution::DiscreteGaussianAlpha {
            alpha,
            mean: None,
            n: None,
        };

        // s_dist = Ternary, m = provided, rough estimation
        let secpar = run_lattice_estimator_cli(
            ring_dim,
            &q,
            &Distribution::Ternary,
            &e_dist,
            Some(m),
            false,
        )?;

        if secpar as u32 >= target_secpar {
            found = Some(found.map_or(mid, |cur| cur.min(mid)));
            // try smaller (more conservative) log_alpha
            hi = mid - 1;
        } else {
            // need larger log_alpha (i.e., larger noise) to increase security
            lo = mid + 1;
        }
    }

    Ok(found.ok_or(SimulatorError::LogAlphaNotFound {
        target_secpar,
        ring_dim: ring_dim.clone(),
        log_q,
        m: m.clone(),
    })?)
}

fn check_correctness(
    target_secpar: u32,
    log_dim: u32,
    crt_bits: u32,
    crt_depth: u32,
    base_bits: u32,
    knapsack_len: u32,
    e_b_log_alpha: i64,
    circuit: &PolyCircuit<DCRTPoly>,
    input_size: usize,
) -> Result<u32, SimulatorError> {
    let ring_dim = BigUint::from(2u32).pow(log_dim);
    let log_q = crt_bits * crt_depth;
    let q = BigUint::from(2u32).pow(log_q);
    let m_g = (crt_bits.div_ceil(base_bits) * crt_depth) as usize;
    let m_b = (m_g + 2) as usize;
    let e_b_sigma = BigDecimal::from_f64(2f64.powf((log_q as i64 - e_b_log_alpha) as f64)).unwrap();
    let secpar_sqrt = BigDecimal::from(target_secpar).sqrt().unwrap();
    let ring_dim_sqrt = BigDecimal::from_biguint(ring_dim.clone(), 0)
        .sqrt()
        .unwrap();
    let base = BigDecimal::from(1 << base_bits);
    let sim_ctx = Arc::new(SimulatorContext::new(
        secpar_sqrt,
        ring_dim_sqrt,
        base,
        m_g as usize,
    ));
    let e_b = PolyMatrixNorm::sample_gauss(sim_ctx.clone(), 1, m_b, e_b_sigma);
    let r_mat = PolyMatrixNorm::new(
        sim_ctx.clone(),
        m_b,
        m_g * input_size,
        BigDecimal::one(),
        Some(m_b - knapsack_len as usize),
    );
    let e_a = &e_b * r_mat;
    let out_wire_norms =
        circuit.simulate_max_h_norm(sim_ctx.clone(), BigDecimal::from(crt_bits), input_size);
    let max_out_wire = out_wire_norms
        .into_iter()
        .max_by(|a, b| a.h_norm.poly_norm.norm.cmp(&b.h_norm.poly_norm.norm))
        .unwrap();
    let (max_h_top, max_h_bottom) = max_out_wire.h_norm.split_rows(m_b);
    let e_after_eval = &e_b * max_h_top + e_a * max_h_bottom;
    let plt_eval = NormPltLweEvaluator::new(sim_ctx.clone(), input_size);
    let e_final = &e_b * plt_eval.preimage1_norm + e_after_eval * plt_eval.preimage2_norm + e_b;
    let q_over_4 = BigDecimal::from_biguint(q, 0) / 4;
    if q_over_4 > e_final.poly_norm.norm {
        Ok(log_dim * m_g as u32)
    } else {
        Err(SimulatorError::NotCorrect {
            e: e_final.poly_norm.norm,
            q_over_4,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mxx::{circuit::PolyCircuit, poly::dcrt::poly::DCRTPoly};

    #[test]
    fn test_bruteforce_params_with_mul() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let ins = circuit.input(2);
        let out_gid = circuit.mul_gate(ins[0], ins[1]);
        circuit.output(vec![out_gid]);

        let params = bruteforce_params(100, 41, (5, 10), (16, 18), (13, 16), circuit, 2);
        assert!(params.is_some());
        println!("params: {:?}", params);
    }
}
