pub mod config;
use crate::config::Config;
use anyhow::Result;
use arithmetic_abe::abe::KeyPolicyABE;
use arithmetic_abe::ciphertext::Ciphertext;
use arithmetic_abe::keys::{FuncSK, MasterPK, MasterSK};
use clap::{Parser, Subcommand};
use keccak_asm::Keccak256;
use mxx::arithmetic::circuit::{ArithGateId, ArithmeticCircuit};
use mxx::{
    matrix::dcrt_poly::DCRTPolyMatrix,
    poly::dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    sampler::{
        PolyTrapdoorSampler, hash::DCRTPolyHashSampler, trapdoor::DCRTPolyTrapdoorSampler,
        uniform::DCRTPolyUniformSampler,
    },
    utils::timed_read,
};
use std::{env, fs};
use std::{path::PathBuf, time::Duration};
use tracing::info;
use tracing_subscriber::{EnvFilter, fmt};

#[derive(Parser, Debug)]
#[command(
    name = "abe",
    version,
    author,
    about = "Key-Policy ABE runner (env-configured)"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Run {
        #[arg(short, long)]
        config: PathBuf,
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
        Commands::Run { config } => run_env_configured(config).await?,
    }

    Ok(())
}

async fn run_env_configured(config: PathBuf) -> Result<()> {
    let contents = fs::read_to_string(&config).unwrap();
    let cfg: Config = toml::from_str(&contents).unwrap();
    let params = DCRTPolyParams::new(
        cfg.ring_dimension,
        cfg.crt_depth,
        cfg.crt_bits,
        cfg.base_bits,
    );

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
    let add_idx = arith.add(ArithGateId::new(0), ArithGateId::new(1)); // a + b
    let mul_idx = arith.mul(add_idx, ArithGateId::new(2)); // (a + b) * c
    let final_idx = arith.sub(mul_idx, ArithGateId::new(0)); // (a + b) * c - a
    arith.output(final_idx);

    // 1) setup
    let (mpk, msk): (
        MasterPK<DCRTPolyMatrix>,
        MasterSK<DCRTPolyMatrix, DCRTPolyTrapdoorSampler>,
    ) = timed_read(
        "setup",
        || abe.setup(params.clone(), cfg.num_inputs, cfg.num_packed_limbs),
        &mut t_setup,
    );
    let config: PathBuf = "keygen".into();
    let dir_path = if config.exists() {
        config
    } else {
        fs::create_dir_all(&config)?;
        config
    };

    info!(target: "abe",  "finished setup");

    // 2) keygen
    let fsk: FuncSK<DCRTPolyMatrix> = abe
        .keygen(
            params.clone(),
            mpk.clone(),
            msk.clone(),
            arith.clone(),
            dir_path,
        )
        .await;

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
