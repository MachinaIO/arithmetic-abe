use mxx::{matrix::PolyMatrix, sampler::PolyTrapdoorSampler};

pub struct MasterPK<M: PolyMatrix> {
    pub num_inputs: usize,
    pub packed_limbs: usize,
    pub seed: [u8; 32],
    b_epsilon: M,
    u: M,
}

impl<M: PolyMatrix> MasterPK<M> {
    pub fn new(num_inputs: usize, packed_limbs: usize, seed: [u8; 32], b_epsilon: M, u: M) -> Self {
        Self {
            num_inputs,
            packed_limbs,
            seed,
            b_epsilon,
            u,
        }
    }
}

pub struct MasterSK<M: PolyMatrix, ST: PolyTrapdoorSampler<M = M>> {
    b_epsilon_trapdoor: ST::Trapdoor,
}

impl<M: PolyMatrix, ST: PolyTrapdoorSampler<M = M>> MasterSK<M, ST> {
    pub fn new(b_epsilon_trapdoor: ST::Trapdoor) -> Self {
        Self { b_epsilon_trapdoor }
    }
}
