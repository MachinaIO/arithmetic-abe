pub mod config;
use crate::config::Config;
use anyhow::Result;
use arithmetic_abe::{
    abe::KeyPolicyABE,
    ciphertext::Ciphertext,
    keys::{FuncSK, MasterPK, MasterSK},
};
use clap::{Parser, Subcommand};
use keccak_asm::Keccak256;
use mxx::{
    arithmetic::circuit::{ArithGateId, ArithmeticCircuit},
    matrix::dcrt_poly::DCRTPolyMatrix,
    poly::dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    sampler::{
        PolyTrapdoorSampler, hash::DCRTPolyHashSampler, trapdoor::DCRTPolyTrapdoorSampler,
        uniform::DCRTPolyUniformSampler,
    },
    utils::timed_read,
};
use std::{env, fs, path::PathBuf, time::Duration};
use tracing::info;
use tracing_subscriber::{EnvFilter, fmt};

#[derive(Parser, Debug)]
#[command(name = "abe", version, author, about = "Key-Policy ABE runner (env-configured)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Run {
        #[arg(short, long)]
        config: PathBuf,

        /// height of the binary tree. Input length should be 2^height
        #[arg(short, long)]
        height: usize,

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
        Commands::Run { config, height, data_dir } => {
            run_env_configured(config, height, data_dir).await?
        }
    }

    Ok(())
}

async fn run_env_configured(config: PathBuf, height: usize, data_dir: PathBuf) -> Result<()> {
    assert_ne!(height, 0);
    let contents = fs::read_to_string(&config).unwrap();
    let cfg: Config = toml::from_str(&contents).unwrap();
    let params =
        DCRTPolyParams::new(cfg.ring_dimension, cfg.crt_depth, cfg.crt_bits, cfg.base_bits);

    let trapdoor_sampler =
        DCRTPolyTrapdoorSampler::new(&params, cfg.trapdoor_sigma.expect("trapdoor sigma exist"));
    let use_packing = false;
    let abe = KeyPolicyABE::<
        DCRTPolyMatrix,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        DCRTPolyUniformSampler,
    >::new(
        cfg.limb_bit_size,
        cfg.crt_bits.div_ceil(cfg.limb_bit_size),
        cfg.crt_depth,
        cfg.num_packed_limbs,
        cfg.d,
        cfg.e_b_sigma,
        use_packing,
        trapdoor_sampler,
    );

    let mut t_setup = Duration::ZERO;
    // let mut t_keygen = Duration::ZERO;
    let mut t_enc = Duration::ZERO;
    let mut t_dec = Duration::ZERO;

    info!(target: "abe",  "starting KeyPolicy ABE");

    let mut arith = ArithmeticCircuit::<DCRTPoly>::setup(
        &params,
        cfg.limb_bit_size,
        cfg.input.len(),
        use_packing,
        true,
    );
    let num_leaves = 1 << (height - 1);
    assert!(
        cfg.input.len() >= num_leaves * 2,
        "Need at least {} inputs for height {} tree",
        num_leaves * 2,
        height
    );

    let mut current_layer = Vec::new();
    for i in 0..num_leaves {
        let left_idx = ArithGateId::new(i * 2);
        let right_idx = ArithGateId::new(i * 2 + 1);
        let mul_gate = arith.mul(left_idx, right_idx);
        current_layer.push(mul_gate);
    }
    for _ in 1..height {
        let mut next_layer = Vec::new();
        let pairs_in_layer = current_layer.len() / 2;

        for i in 0..pairs_in_layer {
            let left = current_layer[i * 2];
            let right = current_layer[i * 2 + 1];
            let mul_gate = arith.mul(left, right);
            next_layer.push(mul_gate);
        }

        current_layer = next_layer;
    }
    assert_eq!(current_layer.len(), 1, "Should have exactly one root gate");
    arith.output(current_layer[0]);

    // 1) setup
    let (mpk, msk): (MasterPK<DCRTPolyMatrix>, MasterSK<DCRTPolyMatrix, DCRTPolyTrapdoorSampler>) =
        timed_read(
            "setup",
            || abe.setup(params.clone(), cfg.num_inputs, cfg.num_packed_limbs),
            &mut t_setup,
        );

    let dir_path = if data_dir.exists() {
        data_dir
    } else {
        fs::create_dir_all(&data_dir)?;
        data_dir
    };

    info!(target: "abe",  "finished setup");

    // 2) keygen
    let fsk: FuncSK<DCRTPolyMatrix> =
        abe.keygen(params.clone(), mpk.clone(), msk.clone(), arith.clone(), dir_path).await;

    // 3) enc
    assert_eq!(cfg.num_inputs, cfg.input.len());
    let ct: Ciphertext<DCRTPolyMatrix> = timed_read(
        "enc",
        || abe.enc(params.clone(), mpk.clone(), &cfg.input, cfg.message),
        &mut t_enc,
    );

    info!(target: "abe", expected_result = cfg.message, "finished encryption");

    // 4) dec
    let bit: bool = timed_read(
        "dec",
        || abe.dec(params.clone(), ct, mpk.clone(), fsk.clone(), arith.clone()),
        &mut t_dec,
    );

    info!(target: "abe", dec_result = bit, "finished decryption");

    Ok(())
}
