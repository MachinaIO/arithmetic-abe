use bigdecimal::{BigDecimal, FromPrimitive, One};
pub use mxx::simulator::lattice_estimator::run_lattice_estimator_cli;
use mxx::{
    arithmetic::circuit::ArithmeticCircuit,
    circuit::PolyCircuit,
    poly::dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    simulator::{
        SimulatorContext,
        lattice_estimator::{Distribution, EstimatorCliError},
        poly_matrix_norm::PolyMatrixNorm,
        wire_norm::NormPltLweEvaluator,
    },
    utils::log_mem,
};
use num_bigint::{BigUint, Sign};
use rayon::{join, prelude::*};
use std::sync::{
    Arc,
    atomic::{AtomicU32, Ordering},
};
use thiserror::Error;
// Logging (replaces println!)
// Configure a logger (e.g., env_logger) in the binary/tests to see output.

#[derive(Debug, Error)]
pub enum SimulatorError {
    /// Error bubbling up from the lattice estimator CLI.
    #[error(transparent)]
    Estimator(#[from] EstimatorCliError),
    /// No knapsack length in [1, max_knapsack_size] reached the target security.
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
    /// No knapsack length in [1, max_knapsack_size] reached the target security.
    #[error(
        "secure knapsack len not found for secpar={target_secpar}, ring_dim={ring_dim}, max_knapsack_size={max_knapsack_size}, q={q}"
    )]
    KnapsackNotFound { target_secpar: u32, ring_dim: BigUint, max_knapsack_size: u32, q: BigUint },
    /// Could not find a log_alpha that reaches the target security.
    #[error(
        "good log_alpha not found for target_secpar={target_secpar}, ring_dim={ring_dim}, log_q={log_q}, m={m}"
    )]
    LogAlphaNotFound { target_secpar: u32, ring_dim: BigUint, log_q: u32, m: BigUint },
    #[error("correctness does not hold: error={e}, q_over_4={q_over_4}")]
    NotCorrect { e: BigDecimal, q_over_4: BigDecimal },
}

// Output (crt_depth, base_bits, log_dim, e_b_sigma, knapsack_size) or None
pub fn bruteforce_params_for_bench_arith_circuit(
    target_secpar: u32,
    crt_bits: u32,
    crt_depth_range: (u32, u32),
    base_bits_range: (u32, u32),
    log_dim_range: (u32, u32),
    num_eval_slots: Option<usize>,
    limb_bit_size: usize,
    height: usize,
    // circuit: PolyCircuit<DCRTPoly>,
) -> Option<(u32, u32, u32, f64, u32)> {
    // (cost, crt_depth, base_bits, log_dim, e_b_sigma, knapsack_size)
    let outputs: Vec<(u32, u32, u32, u32, f64, u32)> =
        (base_bits_range.0..=base_bits_range.1)
            .into_par_iter()
            .flat_map(|base_bits| {
                let mut local = Vec::<(u32, u32, u32, u32, f64, u32)>::new();
                let mut lo = crt_depth_range.0;
                let mut hi = crt_depth_range.1;
                while lo <= hi {
                    let crt_depth = lo + ((hi - lo) / 2);
                    log::info!("base_bits {base_bits} crt_depth {crt_depth}");
                    let (log_dim, e_b_log_alpha, knapsack_size) = match find_min_ring_dim(
                        target_secpar,
                        crt_bits,
                        crt_depth,
                        base_bits,
                        log_dim_range,
                    ) {
                        Ok(result) => result,
                        Err(e) => {
                            log::info!(
                                "Security error with target_secpar = {}, crt_bits = {}, base_bits = {}, crt_depth = {}, limb_bit_size = {}, height = {}: {}",
                                target_secpar, crt_bits, base_bits, crt_depth, limb_bit_size, height, e
                            );
                            // try smaller crt_depth
                            if crt_depth == 0 { break; }
                            hi = crt_depth - 1;
                            continue;
                        }
                    };
                    log::info!(
                        "Found log_dim = {}, e_b_log_alpha = {}, knapsack_size = {}",
                        log_dim,
                        e_b_log_alpha,
                        knapsack_size
                    );
                    let ring_dim = (1 << log_dim) as u32;
                    let params = DCRTPolyParams::new(ring_dim, crt_depth as usize, crt_bits as usize, base_bits);
                    let circuit = ArithmeticCircuit::benchmark_multiplication_tree(&params, limb_bit_size, num_eval_slots.unwrap_or(ring_dim as usize), height,true);
                    log::info!("circuit constructed with crt_depth = {}, log_dim = {}, base_bits = {}, knapsack_size = {}, e_b_log_alpha = {}", crt_depth, log_dim, base_bits, knapsack_size, e_b_log_alpha);
                    log::info!("poly circuit non_free_depth {}",circuit.poly_circuit.non_free_depth());
                    match check_correctness(
                        target_secpar,
                        log_dim,
                        crt_bits,
                        crt_depth,
                        base_bits,
                        knapsack_size,
                        e_b_log_alpha,
                        &circuit.poly_circuit,
                    ) {
                        Ok(cost) => {
                            log::info!(
                                "Found with target_secpar = {}, crt_bits = {}, base_bits = {}, crt_depth = {}, cost = {}",
                                target_secpar, crt_bits, base_bits, crt_depth,  cost
                            );
                            local.push((
                                cost,
                                crt_depth,
                                base_bits,
                                log_dim,
                                2.0f64.powf(
                                    crt_bits as f64 * crt_depth as f64 + e_b_log_alpha as f64,
                                ),
                                knapsack_size,
                            ));
                            // search smaller crt_depth to continue binary search
                            if crt_depth == 0 { break; }
                            hi = crt_depth - 1;
                        }
                        Err(e) => {
                            log::info!(
                                "Correctness error with target_secpar = {}, crt_bits = {}, base_bits = {}, crt_depth = {}: {}",
                                target_secpar, crt_bits, base_bits, crt_depth, e
                            );
                            // try larger crt_depth
                            lo = crt_depth + 1;
                        }
                    }
                }
                local
            })
            .collect();
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
    // Evaluate all candidate log_dim in parallel, then select the minimal feasible one.
    let results: Vec<Result<(u32, i64, u32), SimulatorError>> = (log_dim_range.0..=log_dim_range.1)
        .into_par_iter()
        .map(|log_dim| {
            log::debug!("log_dim {}", log_dim);
            let ring_dim = BigUint::from(2u32).pow(log_dim);
            match check_security(target_secpar, &ring_dim, crt_bits, crt_depth, base_bits) {
                Ok((log_alpha, knapsack_size)) => Ok((log_dim, log_alpha, knapsack_size)),
                Err(e) => Err(e),
            }
        })
        .collect();

    // Pick the smallest log_dim among successes
    if let Some((log_dim, log_alpha, knapsack_size)) =
        results.iter().filter_map(|r| r.as_ref().ok()).min_by(|a, b| a.0.cmp(&b.0)).copied()
    {
        return Ok((log_dim, log_alpha, knapsack_size));
    }

    // If there were no successes, propagate any Estimator error if present
    if let Some(estimator_err) = results.into_iter().find_map(|r| match r {
        Err(SimulatorError::Estimator(e)) => Some(SimulatorError::Estimator(e)),
        _ => None,
    }) {
        return Err(estimator_err);
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
    let log_q = crt_bits * crt_depth;
    let q = BigUint::from(2u32).pow(log_q);
    let m_g = crt_bits.div_ceil(base_bits) * crt_depth;
    let m_b = m_g + 2;
    // The column size of the matrix B (sampled with a trapdoor) is m_b; however, one column is an
    // identity polynomial, so we need to ignore one column. Additionally, one more uniformly
    // random matrix is used for encrypting a message in ABE; thus the total column size for
    // ring-LWE is m_b - 1 + 1 = m_b.
    let (log_alpha_res, knapsack_res) = join(
        || find_log_alpha_for_ring_lwe(target_secpar, ring_dim, log_q, &BigUint::from(m_b)),
        || find_knapsack_size(target_secpar, ring_dim, &q, m_b - 1),
    );
    let log_alpha = log_alpha_res?;
    log::debug!("found log_alpha_res = {log_alpha}");
    let knapsack_size = knapsack_res?;
    log::debug!("found knapsack_size = {knapsack_size}");
    Ok((log_alpha, knapsack_size))
}

/// Returns the smallest `knapsack_size` in [1, max_knapsack_size] whose estimated
/// security is at least `target_secpar`, or an error if estimation fails or none found.
/// - `target_secpar`: required minimum security parameter.
/// - `ring_dim`: base ring dimension.
/// - `q`: modulus (as BigUint).
/// - `max_knapsack_size`: upper bound to search (inclusive).
fn find_knapsack_size(
    target_secpar: u32,
    ring_dim: &BigUint,
    q: &BigUint,
    max_knapsack_size: u32,
) -> Result<u32, SimulatorError> {
    if max_knapsack_size < 2 {
        return Err(SimulatorError::KnapsackNotFound {
            target_secpar,
            ring_dim: ring_dim.clone(),
            max_knapsack_size,
            q: q.clone(),
        });
    }

    let best = AtomicU32::new(0);

    (2..=max_knapsack_size).into_par_iter().try_for_each(
        |knapsack_size| -> Result<(), SimulatorError> {
            let current_best = best.load(Ordering::Relaxed);
            if current_best != 0 && knapsack_size >= current_best {
                return Ok(());
            }

            // Effective LWE dimension n = ring_dim * knapsack_size - ring_dim
            let n = ring_dim * BigUint::from(knapsack_size) - ring_dim;
            // s_dist = Ternary, e_dist = Ternary, m = n, exact = false (rough)
            let secpar = run_lattice_estimator_cli(
                &n,
                q,
                &Distribution::Ternary,
                &Distribution::Ternary,
                Some(&n),
                false,
            )?;
            log::debug!("called estimator {secpar} in find_knapsack_size for {knapsack_size}");

            if secpar as u32 >= target_secpar {
                let mut observed = best.load(Ordering::Acquire);
                while observed == 0 || knapsack_size < observed {
                    match best.compare_exchange(
                        observed,
                        knapsack_size,
                        Ordering::AcqRel,
                        Ordering::Acquire,
                    ) {
                        Ok(_) => break,
                        Err(actual) => {
                            if actual != 0 && knapsack_size >= actual {
                                break;
                            }
                            observed = actual;
                        }
                    }
                }
            }

            Ok(())
        },
    )?;

    let best_value = best.load(Ordering::Relaxed);
    if best_value != 0 {
        Ok(best_value)
    } else {
        Err(SimulatorError::KnapsackNotFound {
            target_secpar,
            ring_dim: ring_dim.clone(),
            max_knapsack_size,
            q: q.clone(),
        })
    }
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
    let mut hi: i64 = 5 - (log_q as i64);
    let mut found: Option<i64> = None;

    while lo <= hi {
        let mid = lo + ((hi - lo) / 2);

        // alpha = sigma/q = 2^{log_alpha}
        let alpha = 2f64.powi(mid as i32); // safe for practical parameter sizes

        let e_dist = Distribution::DiscreteGaussianAlpha { alpha, mean: None, n: None };

        // s_dist = Ternary, m = provided, rough estimation
        let secpar = run_lattice_estimator_cli(
            ring_dim,
            &q,
            &Distribution::Ternary,
            &e_dist,
            Some(m),
            false,
        )?;
        log::debug!("called estimator {secpar} in find_log_alpha_for_ring_lwe");

        if secpar as u32 >= target_secpar {
            found = Some(found.map_or(mid, |cur| cur.min(mid)));
            // try smaller (more conservative) log_alpha
            hi = mid - 1;
        } else {
            // need larger log_alpha (i.e., larger noise) to increase security
            lo = mid + 1;
        }
    }

    found.ok_or(SimulatorError::LogAlphaNotFound {
        target_secpar,
        ring_dim: ring_dim.clone(),
        log_q,
        m: m.clone(),
    })
}

fn check_correctness(
    target_secpar: u32,
    log_dim: u32,
    crt_bits: u32,
    crt_depth: u32,
    base_bits: u32,
    knapsack_size: u32,
    e_b_log_alpha: i64,
    circuit: &PolyCircuit<DCRTPoly>,
) -> Result<u32, SimulatorError> {
    let input_size = circuit.num_input();
    let ring_dim = BigUint::from(2u32).pow(log_dim);
    let log_q = crt_bits * crt_depth;
    let q = BigUint::from(2u32).pow(log_q);
    let m_g = (crt_bits.div_ceil(base_bits) * crt_depth) as usize;
    let m_b = m_g + 2;
    let e_b_sigma = BigDecimal::from_f64(2f64.powf((log_q as i64 - e_b_log_alpha) as f64)).unwrap();
    let secpar_sqrt = BigDecimal::from_u32(target_secpar).unwrap().sqrt().unwrap();
    let ring_dim_sqrt = BigDecimal::from_biguint(ring_dim.clone(), 0).sqrt().unwrap();
    let base = BigDecimal::from_biguint((BigUint::from(1u32)) << base_bits, 0);
    let sim_ctx = Arc::new(SimulatorContext::new(secpar_sqrt, ring_dim_sqrt, base, m_g));
    let e_b = PolyMatrixNorm::sample_gauss(sim_ctx.clone(), 1, m_b, e_b_sigma.clone());

    let r_mat = PolyMatrixNorm::new(
        sim_ctx.clone(),
        m_b,
        m_g * input_size,
        BigDecimal::one(),
        Some(m_b - knapsack_size as usize),
    );
    let e_a = &e_b * &r_mat;
    log::info!("before simulation: e_b = {:?}, e_a = {:?}", e_b, e_a);
    let out_wire_norms = circuit.simulate_max_h_norm(
        sim_ctx.clone(),
        BigDecimal::from_u32(crt_bits).unwrap(),
        input_size,
    );
    log::info!("after simulation");
    let max_out_wire = out_wire_norms
        .into_iter()
        .max_by(|a, b| a.h_norm.poly_norm.norm.cmp(&b.h_norm.poly_norm.norm))
        .unwrap();
    let (max_h_top, max_h_bottom) = max_out_wire.h_norm.split_rows(m_b);
    let e_after_eval = &e_b * max_h_top + e_a * max_h_bottom;
    let plt_eval = NormPltLweEvaluator::new(sim_ctx.clone(), input_size);
    let mut preimage_norm_top = plt_eval.preimage1_norm.clone();
    preimage_norm_top.nrow = m_b;
    preimage_norm_top.ncol = 1;
    let mut preimage_norm_bottom = plt_eval.preimage2_norm.clone();
    preimage_norm_bottom.ncol = 1;
    let e_u = PolyMatrixNorm::sample_gauss(sim_ctx.clone(), 1, 1, e_b_sigma);
    let e_final = &e_b * preimage_norm_top + e_after_eval * preimage_norm_bottom + e_u;
    let q_over_4 = BigDecimal::from_biguint(q, 0) / BigDecimal::from_u32(4).unwrap();
    if q_over_4 > e_final.poly_norm.norm {
        log_mem(format!(
            "q_over_4: {:?}, e_final: {:?}",
            q_over_4.with_scale_round(0, bigdecimal::RoundingMode::Ceiling).to_string(),
            e_final
                .poly_norm
                .norm
                .with_scale_round(0, bigdecimal::RoundingMode::Ceiling)
                .to_string()
        ));
        Ok(log_dim * m_g as u32)
    } else {
        Err(SimulatorError::NotCorrect { e: e_final.poly_norm.norm, q_over_4 })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Initialize logger for test output
    use env_logger;

    #[test]
    fn test_bruteforce_params() {
        // Initialize env_logger once for tests; ignore if already set.
        let _ = env_logger::builder().is_test(true).try_init();
        let params = bruteforce_params_for_bench_arith_circuit(
            100,
            41,
            (2, 4),
            (15, 18),
            (13, 16),
            Some(2),
            2,
            3,
        );
        assert!(params.is_some());
        println!("params: {:?}", params);
    }
}
