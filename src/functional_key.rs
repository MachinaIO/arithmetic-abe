use mxx::matrix::PolyMatrix;
use std::path::PathBuf;

#[derive(Clone)]
pub struct FuncPK {}

#[derive(Clone)]
pub struct FuncSK<M: PolyMatrix> {
    pub a_f: M,
    pub u_f: M,
    pub dir_path: PathBuf,
}
