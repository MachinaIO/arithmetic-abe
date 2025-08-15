use std::str::FromStr;

use num_bigint::BigUint;
use serde::{Deserialize, Deserializer, Serialize, Serializer, de};

fn biguint_to_string<S>(value: &BigUint, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&value.to_str_radix(10))
}

fn biguint_from_string<'de, D>(deserializer: D) -> Result<BigUint, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    BigUint::from_str(&s).map_err(de::Error::custom)
}

fn default_trapdoor_sigma() -> f64 {
    4.578
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub num_inputs: usize,
    pub packed_limbs: usize,
    pub limb_bit_size: usize,
    pub crt_depth: usize,
    pub crt_bits: usize,
    pub d: usize,
    pub p_sigma: f64,
    pub message: u8,
    pub circuit_path: Option<String>,
    #[serde(
        serialize_with = "biguint_to_string",
        deserialize_with = "biguint_from_string"
    )]
    pub switched_modulus: BigUint,
    #[serde(default = "default_trapdoor_sigma")]
    pub trapdoor_sigma: f64,
    /// polynomial ring dimension
    pub ring_dimension: u32,
    /// bit size of the base for the gadget vector and decomposition
    pub base_bits: u32,
    pub input: Vec<u64>,
}
