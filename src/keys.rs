use mxx::{matrix::PolyMatrix, poly::Poly, sampler::PolyTrapdoorSampler};
use std::{
    convert::TryFrom,
    fs::File,
    io::{self, ErrorKind, Read, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

#[derive(Clone)]
pub struct FuncSK<M: PolyMatrix> {
    pub a_f: M,
    pub u_f: M,
    pub dir_path: PathBuf,
}

impl<M: PolyMatrix> FuncSK<M> {
    pub fn new(a_f: M, u_f: M, dir_path: PathBuf) -> Self {
        Self { a_f, u_f, dir_path }
    }

    pub fn write<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut file = File::create(path)?;
        file.write_all(b"FSK0")?;

        let dir_str = self.dir_path.to_string_lossy();
        write_blob(&mut file, dir_str.as_bytes())?;
        write_blob(&mut file, &self.a_f.to_compact_bytes())?;
        write_blob(&mut file, &self.u_f.to_compact_bytes())?;

        Ok(())
    }

    pub fn read<P: AsRef<Path>>(params: &<M::P as Poly>::Params, path: P) -> io::Result<Self> {
        let mut file = File::open(path)?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != b"FSK0" {
            return Err(io::Error::new(ErrorKind::InvalidData, "invalid FuncSK header"));
        }

        let dir_bytes = read_blob(&mut file)?;
        let dir_str = String::from_utf8(dir_bytes)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "dir_path is not valid UTF-8"))?;
        let dir_path = PathBuf::from(dir_str);

        let a_f_bytes = read_blob(&mut file)?;
        let u_f_bytes = read_blob(&mut file)?;

        let a_f = M::from_compact_bytes(params, &a_f_bytes);
        let u_f = M::from_compact_bytes(params, &u_f_bytes);

        Ok(Self { a_f, u_f, dir_path })
    }
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

    pub fn write<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut file = File::create(path)?;
        file.write_all(b"MPK0")?;

        let num_inputs = u64::try_from(self.num_inputs)
            .map_err(|_| io::Error::new(ErrorKind::InvalidInput, "num_inputs exceeds u64"))?;
        file.write_all(&num_inputs.to_le_bytes())?;
        file.write_all(&self.seed)?;

        write_blob(&mut file, &self.b_matrix.to_compact_bytes())?;
        write_blob(&mut file, &self.u.to_compact_bytes())?;

        Ok(())
    }

    pub fn read<P: AsRef<Path>>(params: &<M::P as Poly>::Params, path: P) -> io::Result<Self> {
        let mut file = File::open(path)?;

        let mut magic = [0u8; 4];
        file.read_exact(&mut magic)?;
        if &magic != b"MPK0" {
            return Err(io::Error::new(ErrorKind::InvalidData, "invalid MasterPK header"));
        }

        let num_inputs = read_u64(&mut file)?;
        let num_inputs = usize::try_from(num_inputs)
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "num_inputs too large"))?;

        let mut seed = [0u8; 32];
        file.read_exact(&mut seed)?;

        let b_matrix_bytes = read_blob(&mut file)?;
        let u_bytes = read_blob(&mut file)?;

        let b_matrix = Arc::new(M::from_compact_bytes(params, &b_matrix_bytes));
        let u = M::from_compact_bytes(params, &u_bytes);

        Ok(Self { num_inputs, seed, b_matrix, u })
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

fn write_blob<W: Write>(writer: &mut W, bytes: &[u8]) -> io::Result<()> {
    let len = u64::try_from(bytes.len())
        .map_err(|_| io::Error::new(ErrorKind::InvalidInput, "blob length exceeds u64"))?;
    writer.write_all(&len.to_le_bytes())?;
    writer.write_all(bytes)?;
    Ok(())
}

fn read_blob<R: Read>(reader: &mut R) -> io::Result<Vec<u8>> {
    let len = read_u64(reader)? as usize;
    let mut buf = vec![0u8; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn read_u64<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}
