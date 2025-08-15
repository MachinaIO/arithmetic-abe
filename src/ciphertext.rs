use mxx::{bgg::encoding::BggEncoding, matrix::PolyMatrix};

pub struct Ciphertext<M: PolyMatrix> {
    pub bgg_encodings: Vec<BggEncoding<M>>,
    pub c_b_epsilon: M,
    pub c_u: M::P,
}
