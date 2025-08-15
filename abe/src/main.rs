pub mod config;
use crate::config::Config;
use anyhow::Result;
use arithmetic_abe::abe::KeyPolicyABE;
use arithmetic_abe::ciphertext::Ciphertext;
use arithmetic_abe::circuit::ArithmeticCircuit;
use arithmetic_abe::keys::{FuncSK, MasterPK, MasterSK};
use clap::{Parser, Subcommand};
use keccak_asm::Keccak256;
use mxx::{
    circuit::PolyCircuit,
    matrix::dcrt_poly::DCRTPolyMatrix,
    poly::dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    sampler::{
        PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler, hash::DCRTPolyHashSampler,
        trapdoor::DCRTPolyTrapdoorSampler, uniform::DCRTPolyUniformSampler,
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

fn main() -> Result<()> {
    init_logging();
    let cli = Cli::parse();

    match cli.command {
        Commands::Run { config } => run_env_configured(config)?,
    }

    Ok(())
}

fn run_env_configured(config: PathBuf) -> Result<()> {
    let contents = fs::read_to_string(&config).unwrap();
    let cfg: Config = toml::from_str(&contents).unwrap();
    let params = DCRTPolyParams::new(
        cfg.ring_dimension,
        cfg.crt_depth,
        cfg.crt_bits,
        cfg.base_bits,
    );
    let base_circuit: ArithmeticCircuit<DCRTPoly> = ArithmeticCircuit {
        limb_bit_size: cfg.limb_bit_size,
        crt_bits: cfg.crt_bits,
        crt_depth: cfg.crt_depth,
        num_crt_limbs: cfg.num_inputs,
        packed_limbs: cfg.packed_limbs,
        num_inputs: cfg.num_inputs,
        original_circuit: PolyCircuit::new(),
    };

    let uniform_sampler = DCRTPolyUniformSampler::new();
    let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new();
    let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, cfg.trapdoor_sigma);
    let abe = KeyPolicyABE::<
        DCRTPolyMatrix,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        DCRTPolyUniformSampler,
    > {
        limb_bit_size: cfg.limb_bit_size,
        num_crt_limbs: cfg.crt_bits.div_ceil(cfg.limb_bit_size),
        crt_depth: cfg.crt_depth,
        packed_limb: cfg.packed_limbs,
        d: cfg.d,
        hash_sampler,
        trapdoor_sampler,
        uniform_sampler,
        p_sigma: cfg.p_sigma,
    };

    let mut t_setup = Duration::ZERO;
    let mut t_keygen = Duration::ZERO;
    let mut t_enc = Duration::ZERO;
    let mut t_dec = Duration::ZERO;

    info!(target: "abe",  "starting KeyPolicy ABE");

    let arith = base_circuit.clone();

    // 1) setup
    let (mpk, msk): (
        MasterPK<DCRTPolyMatrix>,
        MasterSK<DCRTPolyMatrix, DCRTPolyTrapdoorSampler>,
    ) = timed_read(
        "setup",
        || abe.setup(params.clone(), cfg.num_inputs, cfg.packed_limbs),
        &mut t_setup,
    );

    // 3) enc
    assert_eq!(cfg.num_inputs, cfg.input.len());
    // let inputs = make_inputs::<DCRTPoly>(&params, cfg.input);
    let msg_bit = cfg.message != 0;
    let ct: Ciphertext<DCRTPolyMatrix> = timed_read(
        "enc",
        || abe.enc(params.clone(), mpk.clone(), &cfg.input, msg_bit),
        &mut t_enc,
    );

    // 2) keygen
    let fsk: FuncSK<DCRTPolyMatrix> = timed_read(
        "keygen",
        || abe.keygen(params.clone(), mpk.clone(), msk.clone(), arith.clone()),
        &mut t_keygen,
    );

    // 4) dec
    let bit: bool = timed_read(
        "dec",
        || abe.dec(params.clone(), ct, mpk.clone(), fsk.clone(), arith.clone()),
        &mut t_dec,
    );

    info!(target: "abe", dec_result = bit, "decryption finished");

    Ok(())
}
