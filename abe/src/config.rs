use serde::{Deserialize, Serialize};

fn default_trapdoor_sigma() -> Option<f64> {
    Some(4.578)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunConfig {
    pub config_id: String,
    pub target_secpar: u32,
    pub crt_depth: u32,
    pub crt_bits: u32,
    pub ring_dimension: u32,
    pub num_eval_slots: Option<usize>,
    pub knapsack_size: Option<usize>,
    pub e_b_sigma: f64,
    #[serde(default = "default_trapdoor_sigma")]
    pub trapdoor_sigma: Option<f64>,
    /// bit size of the base for the gadget vector and decomposition
    pub base_bits: u32,
    pub l1_moduli_bits: usize,
    pub scale: u64,
    pub arith_input_size: usize,
    pub arith_height: u32,
    // #[serde(
    //     serialize_with = "matrix_biguint_to_string",
    //     deserialize_with = "matrix_biguint_from_string"
    // )]
    // pub input_rand_seed: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    pub target_secpar: u32,
    pub crt_bits: u32,
    pub crt_depth_min: u32,
    pub crt_depth_max: u32,
    pub base_bits_min: u32,
    pub base_bits_max: u32,
    pub log_dim_min: u32,
    pub log_dim_max: u32,
    pub num_eval_slots: Option<usize>,
    pub l1_moduli_bits: usize,
    pub scale: u64,
    pub height: usize,
}
