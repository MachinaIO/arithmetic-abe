pub use mxx::simulator::lattice_estimator::run_lattice_estimator_cli;

use mxx::simulator::lattice_estimator::{Distribution, EstimatorCliError};
use num_bigint::BigUint;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SimulatorError {
    /// Error bubbling up from the lattice estimator CLI.
    #[error(transparent)]
    Estimator(#[from] EstimatorCliError),
    /// No knapsack length in [1, max_knapsack_len] reached the target security.
    #[error(
        "secure knapsack len not found for secpar={target_secpar}, ring_dim={ring_dim}, max_knapsack_len={max_knapsack_len}, q={q}"
    )]
    /// Could not find a knapsack_len that reaches the target security.
    KnapsackNotFound {
        target_secpar: u64,
        ring_dim: BigUint,
        max_knapsack_len: usize,
        q: BigUint,
    },
    /// Could not find a log_alpha that reaches the target security.
    #[error(
        "good log_alpha not found for target_secpar={target_secpar}, ring_dim={ring_dim}, log_q={log_q}, m={m}"
    )]
    LogAlphaNotFound {
        target_secpar: u64,
        ring_dim: BigUint,
        log_q: u64,
        m: BigUint,
    },
}

fn check_ring_dim(target_secpar: u64, ring_dim: &BigUint, crt_bits: usize, crt_depth: usize) {}

/// Returns the smallest `knapsack_len` in [1, max_knapsack_len] whose estimated
/// security is at least `target_secpar`, or an error if estimation fails or none found.
/// - `target_secpar`: required minimum security parameter.
/// - `ring_dim`: base ring dimension.
/// - `q`: modulus (as BigUint).
/// - `max_knapsack_len`: upper bound to search (inclusive).
fn find_knapsack_len(
    target_secpar: u64,
    ring_dim: &BigUint,
    q: &BigUint,
    max_knapsack_len: usize,
) -> Result<usize, SimulatorError> {
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

        if secpar >= target_secpar {
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
    target_secpar: u64,
    ring_dim: &BigUint,
    log_q: u64,
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

        if secpar >= target_secpar {
            found = Some(found.map_or(mid, |cur| cur.min(mid)));
            // try smaller (more conservative) log_alpha
            hi = mid - 1;
        } else {
            // need larger log_alpha (i.e., larger noise) to increase security
            lo = mid + 1;
        }
    }

    if let Some(ans) = found {
        Ok(ans)
    } else {
        Err(SimulatorError::LogAlphaNotFound {
            target_secpar,
            ring_dim: ring_dim.clone(),
            log_q,
            m: m.clone(),
        })
    }
}
