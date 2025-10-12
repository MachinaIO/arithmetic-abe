pub mod config;
use crate::config::{RunConfig, SimConfig};
use anyhow::{Context, Result, ensure};
use arithmetic_abe::{
    abe::KeyPolicyABE,
    ciphertext::Ciphertext,
    keys::{FuncSK, MasterPK, MasterSK},
    simulator::bruteforce_params_for_bench_arith_circuit,
};
use chrono::Local;
use clap::{Parser, Subcommand};
use keccak_asm::Keccak256;
use mxx::{
    arithmetic::circuit::ArithmeticCircuit,
    matrix::dcrt_poly::DCRTPolyMatrix,
    poly::{
        PolyParams,
        dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    },
    sampler::{
        PolyTrapdoorSampler, hash::DCRTPolyHashSampler, trapdoor::DCRTPolyTrapdoorSampler,
        uniform::DCRTPolyUniformSampler,
    },
    utils::{log_mem, timed_read, timed_read_async},
};
use num_bigint::BigUint;
use std::{env, fs, path::PathBuf, time::Duration};
use tracing_subscriber::{EnvFilter, fmt};

#[derive(Parser, Debug)]
#[command(name = "abe", version, author, about = "Key-Policy ABE runner (env-configured)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    BenchSim {
        #[arg(
            short,
            long,
            value_name = "PATH",
            help = "Path to a TOML file containing SimConfig parameters"
        )]
        config: PathBuf,
    },
    BenchRunOffline {
        #[arg(short, long)]
        config: PathBuf,

        #[arg(short, long)]
        data_dir: PathBuf,
    },
}

fn init_logging() {
    let level = env::var("RUST_LOG").ok().unwrap_or("info".into());
    let env_filter = EnvFilter::try_new(level).unwrap_or_else(|_| EnvFilter::new("info"));
    fmt().with_env_filter(env_filter).with_target(true).init();
}

#[tokio::main]
async fn main() -> Result<()> {
    init_logging();
    let cli = Cli::parse();

    match cli.command {
        Commands::BenchSim { config } => {
            log_mem(format!("Loading simulator config: path={}", config.display()));
            let config_prefix = config
                .file_name()
                .map(|os| os.to_string_lossy().into_owned())
                .unwrap_or_else(|| "sim-config".to_string());
            let contents = fs::read_to_string(&config).with_context(|| {
                format!("failed to read simulator config from {}", config.display())
            })?;
            let sim_config: SimConfig = toml::from_str(&contents)
                .with_context(|| format!("failed to parse SimConfig from {}", config.display()))?;
            run_bench_sim(sim_config, config_prefix)?;
        }
        Commands::BenchRunOffline { config, data_dir } => {
            log_mem(format!("Loading run config: path={}", config.display()));
            let contents = fs::read_to_string(&config)
                .with_context(|| format!("failed to read run config from {}", config.display()))?;
            let run_config: RunConfig = toml::from_str(&contents)
                .with_context(|| format!("failed to parse RunConfig from {}", config.display()))?;
            run_bench_offline(run_config, data_dir).await?;
        }
    }

    Ok(())
}

fn run_bench_sim(config: SimConfig, config_prefix: String) -> Result<()> {
    ensure!(config.crt_depth_min <= config.crt_depth_max, "invalid CRT depth range: min > max");
    ensure!(config.base_bits_min <= config.base_bits_max, "invalid base bits range: min > max");
    ensure!(config.log_dim_min <= config.log_dim_max, "invalid log dim range: min > max");

    let SimConfig {
        target_secpar,
        crt_bits,
        crt_depth_min,
        crt_depth_max,
        base_bits_min,
        base_bits_max,
        log_dim_min,
        log_dim_max,
        num_eval_slots,
        limb_bit_size,
        height,
    } = config;

    log_mem(format!(
        "Starting benchmark parameter search: target_secpar={}, crt_bits={}, crt_depth_range=({}-{}), base_bits_range=({}-{}), log_dim_range=({}-{}), num_eval_slots={:?}, limb_bit_size={}, height={}, config_prefix={}",
        target_secpar,
        crt_bits,
        crt_depth_min,
        crt_depth_max,
        base_bits_min,
        base_bits_max,
        log_dim_min,
        log_dim_max,
        num_eval_slots,
        limb_bit_size,
        height,
        config_prefix
    ));

    let params = bruteforce_params_for_bench_arith_circuit(
        target_secpar,
        crt_bits,
        (crt_depth_min, crt_depth_max),
        (base_bits_min, base_bits_max),
        (log_dim_min, log_dim_max),
        config.num_eval_slots,
        limb_bit_size,
        height,
    )
    .context("unable to find parameters for benchmark arithmetic circuit")?;

    let (crt_depth, base_bits, log_dim, e_b_sigma, knapsack_size) = params;
    let ring_dimension =
        1u32.checked_shl(log_dim).context("log_dim is too large for u32 ring dimension")?;
    let knapsack_size = knapsack_size as usize;
    let arith_height = height as u32;
    let arith_input_size = 1usize
        .checked_shl(arith_height)
        .context("arith_height is too large for usize input size")?;

    log_mem(format!(
        "Benchmark parameter search succeeded: crt_depth={}, base_bits={}, log_dim={}, e_b_sigma={}, knapsack_size={}, ring_dimension={}, arith_input_size={}",
        crt_depth, base_bits, log_dim, e_b_sigma, knapsack_size, ring_dimension, arith_input_size
    ));

    let config_id = format!("{}_{}", config_prefix, Local::now().format("%Y%m%d-%H%M%S"));
    let run_config = RunConfig {
        config_id: config_id.clone(),
        target_secpar,
        crt_depth,
        crt_bits,
        ring_dimension,
        knapsack_size: Some(knapsack_size),
        e_b_sigma,
        trapdoor_sigma: Some(4.578),
        base_bits,
        num_eval_slots: config.num_eval_slots,
        limb_bit_size,
        arith_input_size,
        arith_height,
    };

    let params_dir = PathBuf::from("abe").join("run_configs");
    fs::create_dir_all(&params_dir).with_context(|| {
        format!("failed to create params directory at {}", params_dir.display())
    })?;
    let output_path = params_dir.join(format!("{}.params.toml", config_id));
    let toml =
        toml::to_string_pretty(&run_config).context("failed to serialize RunConfig into TOML")?;
    fs::write(&output_path, toml)
        .with_context(|| format!("failed to write config file to {}", output_path.display()))?;

    log_mem(format!("Wrote benchmark config: path={}", output_path.display()));

    Ok(())
}

async fn run_bench_offline(config: RunConfig, data_dir: PathBuf) -> Result<()> {
    let params = DCRTPolyParams::new(
        config.ring_dimension,
        config.crt_depth as usize,
        config.crt_bits as usize,
        config.base_bits,
    );
    let trapdoor_sampler =
        DCRTPolyTrapdoorSampler::new(&params, config.trapdoor_sigma.expect("trapdoor sigma exist"));
    let abe = KeyPolicyABE::<
        DCRTPolyMatrix,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        DCRTPolyUniformSampler,
    >::new(
        config.limb_bit_size,
        &params,
        config.num_eval_slots,
        config.knapsack_size,
        config.e_b_sigma,
        trapdoor_sampler,
    );
    let mut t_setup = Duration::ZERO;
    let mut t_keygen = Duration::ZERO;

    log_mem("starting KeyPolicy ABE");
    log_mem("start building arithmetic circuit");
    let arith_circuit = ArithmeticCircuit::<DCRTPoly>::benchmark_multiplication_tree(
        &params,
        config.limb_bit_size,
        config.ring_dimension as usize,
        config.arith_height as usize,
        false,
    );
    log_mem("finished building arithmetic circuit");

    // 1) setup
    log_mem("starting setup");
    let (mpk, msk): (MasterPK<DCRTPolyMatrix>, MasterSK<DCRTPolyMatrix, DCRTPolyTrapdoorSampler>) =
        timed_read("setup", || abe.setup(params.clone(), config.arith_input_size), &mut t_setup);
    log_mem("finished setup");

    let dir_path = if data_dir.exists() {
        data_dir
    } else {
        fs::create_dir_all(&data_dir)?;
        data_dir
    };
    // 2) keygen
    log_mem("starting keygen");
    let fsk: FuncSK<DCRTPolyMatrix> = timed_read_async(
        "keygen",
        || abe.keygen(params.clone(), mpk.clone(), msk.clone(), arith_circuit, dir_path.clone()),
        &mut t_keygen,
    )
    .await;
    log_mem("finished keygen");

    log_mem("starting writing mpk and fsk files");
    mpk.write(dir_path.join(format!("{}.mpk", config.config_id)))?;
    fsk.write(dir_path.join(format!("{}.fsk", config.config_id)))?;
    log_mem("finished writing mpk and fsk files");
    Ok(())
}

async fn run_bench_online(config: RunConfig, data_dir: PathBuf) -> Result<()> {
    let params = DCRTPolyParams::new(
        config.ring_dimension,
        config.crt_depth as usize,
        config.crt_bits as usize,
        config.base_bits,
    );
    let trapdoor_sampler =
        DCRTPolyTrapdoorSampler::new(&params, config.trapdoor_sigma.expect("trapdoor sigma exist"));
    let abe = KeyPolicyABE::<
        DCRTPolyMatrix,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        DCRTPolyUniformSampler,
    >::new(
        config.limb_bit_size,
        &params,
        config.num_eval_slots,
        config.knapsack_size,
        config.e_b_sigma,
        trapdoor_sampler,
    );
    let mut t_read_mpk = Duration::ZERO;
    let mut t_enc = Duration::ZERO;
    let mut t_read_fsk = Duration::ZERO;
    let mut t_dec = Duration::ZERO;

    log_mem("starting KeyPolicy ABE");
    log_mem("start building arithmetic circuit");
    let arith_circuit = ArithmeticCircuit::<DCRTPoly>::benchmark_multiplication_tree(
        &params,
        config.limb_bit_size,
        config.ring_dimension as usize,
        config.arith_height as usize,
        false,
    );
    log_mem("finished building arithmetic circuit");

    let ring_dim = params.ring_dimension() as usize;
    // 3) enc
    log_mem("starting enc");
    let mpk = timed_read(
        "read mpk",
        || {
            MasterPK::<DCRTPolyMatrix>::read(
                &params,
                data_dir.join(format!("{}.mpk", config.config_id)),
            )
            .expect("failed to read mpk file")
        },
        &mut t_read_mpk,
    );
    let ct: Ciphertext<DCRTPolyMatrix> = timed_read(
        "enc",
        || {
            abe.enc(
                params.clone(),
                mpk,
                &vec![vec![BigUint::ZERO; ring_dim]; config.arith_input_size],
                &vec![true; ring_dim],
            )
        },
        &mut t_enc,
    );
    log_mem("finished enc");
    // 4) dec
    log_mem("starting dec");
    t_read_mpk = Duration::ZERO;
    let mpk = timed_read(
        "read mpk",
        || {
            MasterPK::<DCRTPolyMatrix>::read(
                &params,
                data_dir.join(format!("{}.mpk", config.config_id)),
            )
            .expect("failed to read mpk file")
        },
        &mut t_read_mpk,
    );
    let fsk = timed_read(
        "read fsk",
        || {
            FuncSK::<DCRTPolyMatrix>::read(
                &params,
                data_dir.join(format!("{}.fsk", config.config_id)),
            )
            .expect("failed to read fsk file")
        },
        &mut t_read_fsk,
    );
    let bit: bool =
        timed_read("dec", || abe.dec(params.clone(), ct, mpk, fsk, arith_circuit), &mut t_dec);
    log_mem(format!("finished decryption: result={}", bit));
    Ok(())
}

// async fn run_env_configured(config: PathBuf, height: usize, data_dir: PathBuf) -> Result<()> {
//     assert_ne!(height, 0);
//     let contents = fs::read_to_string(&config).unwrap();
//     let cfg: RunConfig = toml::from_str(&contents).unwrap();
//     let params =
//         DCRTPolyParams::new(cfg.ring_dimension, cfg.crt_depth, cfg.crt_bits, cfg.base_bits);

//     let trapdoor_sampler =
//         DCRTPolyTrapdoorSampler::new(&params, cfg.trapdoor_sigma.expect("trapdoor sigma exist"));
//     let use_packing = false;
//     let abe = KeyPolicyABE::<
//         DCRTPolyMatrix,
//         DCRTPolyHashSampler<Keccak256>,
//         DCRTPolyTrapdoorSampler,
//         DCRTPolyUniformSampler,
//     >::new(
//         cfg.limb_bit_size,
//         cfg.crt_bits.div_ceil(cfg.limb_bit_size),
//         cfg.crt_depth,
//         cfg.knapsack_size,
//         cfg.e_b_sigma,
//         use_packing,
//         trapdoor_sampler,
//     );

//     let mut t_setup = Duration::ZERO;
//     // let mut t_keygen = Duration::ZERO;
//     let mut t_enc = Duration::ZERO;
//     let mut t_dec = Duration::ZERO;

//     info!(target: "abe",  "starting KeyPolicy ABE");

//     let num_inputs = cfg.input.len();
//     let mut arith = ArithmeticCircuit::<DCRTPoly>::setup(
//         &params,
//         cfg.limb_bit_size,
//         cfg.ring_dimension as usize,
//         num_inputs,
//         use_packing,
//         true,
//     );
//     let num_leaves = 1 << (height - 1);
//     assert!(
//         cfg.input.len() >= num_leaves * 2,
//         "Need at least {} inputs for height {} tree",
//         num_leaves * 2,
//         height
//     );
//     info!("setup done");

//     let mut current_layer = Vec::new();
//     for i in 0..num_leaves {
//         let left_idx = ArithGateId::new(i * 2);
//         let right_idx = ArithGateId::new(i * 2 + 1);
//         let mul_gate = arith.mul(left_idx, right_idx);
//         current_layer.push(mul_gate);
//     }
//     for _ in 1..height {
//         let mut next_layer = Vec::new();
//         let pairs_in_layer = current_layer.len() / 2;

//         for i in 0..pairs_in_layer {
//             let left = current_layer[i * 2];
//             let right = current_layer[i * 2 + 1];
//             let mul_gate = arith.mul(left, right);
//             next_layer.push(mul_gate);
//         }

//         current_layer = next_layer;
//     }
//     assert_eq!(current_layer.len(), 1, "Should have exactly one root gate");
//     arith.output(current_layer[0]);

//     // 1) setup
//     let (mpk, msk): (MasterPK<DCRTPolyMatrix>, MasterSK<DCRTPolyMatrix, DCRTPolyTrapdoorSampler>)
// =         timed_read("setup", || abe.setup(params.clone(), num_inputs), &mut t_setup);

//     let dir_path = if data_dir.exists() {
//         data_dir
//     } else {
//         fs::create_dir_all(&data_dir)?;
//         data_dir
//     };

//     info!(target: "abe",  "finished setup");

//     // 2) keygen
//     let fsk: FuncSK<DCRTPolyMatrix> =
//         abe.keygen(params.clone(), mpk.clone(), msk.clone(), arith.clone(), dir_path).await;

//     // 3) enc
//     let ct: Ciphertext<DCRTPolyMatrix> = timed_read(
//         "enc",
//         || abe.enc(params.clone(), mpk.clone(), &cfg.input, &cfg.message),
//         &mut t_enc,
//     );

//     // 4) dec
//     let bit: bool = timed_read(
//         "dec",
//         || abe.dec(params.clone(), ct, mpk.clone(), fsk.clone(), arith.clone()),
//         &mut t_dec,
//     );

//     info!(target: "abe", dec_result = bit, "finished decryption");

//     Ok(())
// }
