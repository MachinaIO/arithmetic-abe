use std::str::FromStr;

use num_bigint::BigUint;
use serde::{Deserialize, Deserializer, Serialize, Serializer, de};

fn matrix_biguint_to_string<S>(values: &Vec<Vec<BigUint>>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let strings: Vec<Vec<String>> =
        values.iter().map(|row| row.iter().map(|v| v.to_str_radix(10)).collect()).collect();
    strings.serialize(serializer)
}

fn matrix_biguint_from_string<'de, D>(deserializer: D) -> Result<Vec<Vec<BigUint>>, D::Error>
where
    D: Deserializer<'de>,
{
    let strings: Vec<Vec<String>> = Vec::deserialize(deserializer)?;
    strings
        .into_iter()
        .map(|row| {
            row.into_iter().map(|s| BigUint::from_str(&s).map_err(de::Error::custom)).collect()
        })
        .collect()
}

fn default_trapdoor_sigma() -> Option<f64> {
    Some(4.578)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub limb_bit_size: usize,
    pub crt_depth: usize,
    pub crt_bits: usize,
    pub knapsack_size: Option<usize>,
    pub e_b_sigma: f64,
    pub message: Vec<bool>,
    #[serde(default = "default_trapdoor_sigma")]
    pub trapdoor_sigma: Option<f64>,
    /// polynomial ring dimension
    pub ring_dimension: u32,
    /// bit size of the base for the gadget vector and decomposition
    pub base_bits: u32,
    #[serde(
        serialize_with = "matrix_biguint_to_string",
        deserialize_with = "matrix_biguint_from_string"
    )]
    pub input: Vec<Vec<BigUint>>,
}
