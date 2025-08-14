use std::path::PathBuf;

use crate::circuit::ArithmeticCircuit;
use mxx::matrix::PolyMatrix;

pub struct FuncPK {}

pub struct FuncSK<M: PolyMatrix> {
    pub arith_circuit: ArithmeticCircuit<M::P>,
    pub a_f: M,
    pub u_f: M,
    pub dir_path: PathBuf,
}
