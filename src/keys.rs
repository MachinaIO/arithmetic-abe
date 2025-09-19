use mxx::{matrix::PolyMatrix, sampler::PolyTrapdoorSampler};
use std::{path::PathBuf, sync::Arc};

#[derive(Clone)]
pub struct FuncSK<M: PolyMatrix> {
    pub a_f: M,
    pub u_f: M,
    pub dir_path: PathBuf,
}

#[derive(Clone)]
pub struct MasterPK<M: PolyMatrix> {
    pub num_inputs: usize,
    pub seed: [u8; 32],
    pub b_matrix: Arc<M>,
    pub u: M,
}

impl<M: PolyMatrix> MasterPK<M> {
    pub fn new(num_inputs: usize, seed: [u8; 32], b_matrix: Arc<M>, u: M) -> Self {
        Self { num_inputs, seed, b_matrix, u }
    }
}

#[derive(Clone)]
pub struct MasterSK<M: PolyMatrix, ST: PolyTrapdoorSampler<M = M>> {
    pub b_trapdoor: Arc<ST::Trapdoor>,
}

impl<M: PolyMatrix, ST: PolyTrapdoorSampler<M = M>> MasterSK<M, ST> {
    pub fn new(b_trapdoor: Arc<ST::Trapdoor>) -> Self {
        Self { b_trapdoor }
    }
}
