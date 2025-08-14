use crate::circuit::ArithmeticCircuit;
use mxx::matrix::PolyMatrix;

pub struct FuncPK {}

pub struct FuncSK<M: PolyMatrix> {
    pub arith_circuit: ArithmeticCircuit<M::P>,
}
