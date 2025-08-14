use std::path::PathBuf;

use crate::circuit::ArithmeticCircuit;
use mxx::matrix::PolyMatrix;

#[derive(Clone)]
pub struct FuncPK {}

#[derive(Clone)]
pub struct FuncSK<M: PolyMatrix> {
    pub arith_circuit: ArithmeticCircuit<M::P>,
    pub a_f: M,
    pub u_f: M,
    pub dir_path: PathBuf,
}
