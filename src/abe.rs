use crate::{
    ciphertext::Ciphertext,
    keys::FuncSK,
    keys::{MasterPK, MasterSK},
};
use mxx::element::PolyElem;
use mxx::gadgets::packed_crt::biguints_to_packed_crt_polys;
use mxx::{
    arithmetic::circuit::ArithmeticCircuit,
    gadgets::{crt::num_limbs_of_crt_poly, packed_crt::num_packed_crt_poly},
};
use mxx::{
    bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
};
use mxx::{gadgets::crt::biguint_to_crt_poly, lookup::lwe_eval::LweBggEncodingPltEvaluator};
use num_bigint::BigUint;
use std::path::PathBuf;
use std::sync::Arc;

const TAG_BGG_PUBKEY: &[u8] = b"BGG_PUBKEY";

pub struct KeyPolicyABE<
    M: PolyMatrix + 'static,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
> {
    pub p_sigma: f64,
    pub limb_bit_size: usize,
    pub num_crt_limbs: usize,
    pub crt_depth: usize,
    pub packed_limb: usize,
    pub d: usize,
    pub use_packing: bool,
    pub hash_sampler: SH,
    pub trapdoor_sampler: ST,
    pub uniform_sampler: SU,
}

impl<
    M: PolyMatrix + 'static,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
> KeyPolicyABE<M, SH, ST, SU>
{
    pub fn setup(
        &self,
        params: <M::P as Poly>::Params,
        num_inputs: usize,
        packed_limbs: usize,
    ) -> (MasterPK<M>, MasterSK<M, ST>) {
        let seed: [u8; 32] = rand::random();
        let (b_epsilon_trapdoor, b_epsilon) = self.trapdoor_sampler.trapdoor(&params, self.d + 1);
        let b_epsilon_trapdoor = Arc::new(b_epsilon_trapdoor);
        let b_epsilon = Arc::new(b_epsilon);
        let u = self
            .uniform_sampler
            .sample_uniform(&params, self.d, 1, DistType::BitDist);
        let mpk = MasterPK::new(num_inputs, packed_limbs, seed, b_epsilon, u);
        let msk = MasterSK::new(b_epsilon_trapdoor);
        (mpk, msk)
    }

    pub fn enc(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        inputs: &[BigUint],
        message: bool,
    ) -> Ciphertext<M> {
        let num_inputs = inputs.len();
        let s = &self
            .uniform_sampler
            .sample_uniform(&params, 1, self.d, DistType::BitDist);
        let uniform_sampler = SU::new();
        let bgg_encoding_sampler =
            BGGEncodingSampler::new(&params, &s.get_row(0), uniform_sampler, self.p_sigma);
        let plaintexts = if self.use_packing {
            biguints_to_packed_crt_polys(self.limb_bit_size, &params, inputs)
        } else {
            inputs
                .into_iter()
                .flat_map(|input| biguint_to_crt_poly(self.limb_bit_size, &params, input))
                .collect()
        };
        let num_given_input_polys = if self.use_packing {
            num_packed_crt_poly::<M::P>(self.limb_bit_size, &params, num_inputs)
        } else {
            num_inputs * num_limbs_of_crt_poly::<M::P>(self.limb_bit_size, &params)
        };
        let reveal_plaintexts = vec![true; num_given_input_polys + 1];
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let pubkeys = bgg_pubkey_sampler.sample(&params, TAG_BGG_PUBKEY, &reveal_plaintexts);
        let bgg_encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let e_cu = &self.uniform_sampler.sample_poly(
            &params,
            &DistType::GaussDist {
                sigma: self.p_sigma,
            },
        );
        let c_b_epsilon_error = &self.uniform_sampler.sample_uniform(
            &params,
            1,
            self.d * (2 + params.modulus_digits()),
            DistType::GaussDist {
                sigma: self.p_sigma,
            },
        );
        let c_b_epsilon = s.clone() * mpk.b_epsilon.as_ref() + c_b_epsilon_error;
        let scale = M::P::from_elem_to_constant(
            &params,
            &(<M::P as Poly>::Elem::half_q(&params.modulus())
                * <M::P as Poly>::Elem::new(message, params.modulus())),
        );
        let c_u = (s.clone() * mpk.u.clone()).get_row(0)[0].clone() + e_cu + scale;

        Ciphertext {
            bgg_encodings,
            c_b_epsilon,
            c_u,
        }
    }

    pub async fn keygen(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, ST>,
        arith_circuit: ArithmeticCircuit<M::P>,
    ) -> FuncSK<M> {
        let dir_path: PathBuf = "keygen".into();
        let result = arith_circuit
            .evaluate_with_bgg_pubkey::<M, SH, ST, SU>(
                &params,
                mpk.seed,
                dir_path.clone(),
                self.d,
                mpk.b_epsilon.clone(),
                msk.b_epsilon_trapdoor.clone(),
                self.trapdoor_sampler.clone(),
            )
            .await;
        let a_f = result[0].clone().matrix;
        let u_f = self.trapdoor_sampler.preimage_extend(
            &params,
            &msk.b_epsilon_trapdoor,
            &mpk.b_epsilon,
            &a_f,
            &mpk.u,
        );
        assert_eq!(result.len(), 1);
        FuncSK { a_f, u_f, dir_path }
    }

    pub fn dec(
        &self,
        params: <M::P as Poly>::Params,
        ct: Ciphertext<M>,
        mpk: MasterPK<M>,
        fsk: FuncSK<M>,
        arith_circuit: ArithmeticCircuit<M::P>,
    ) -> bool {
        let encodings = &ct.bgg_encodings[..];
        let dir_path: PathBuf = fsk.dir_path;
        let bgg_evaluator =
            LweBggEncodingPltEvaluator::<M, SH>::new(mpk.seed, dir_path, ct.c_b_epsilon.clone());
        let result = arith_circuit.poly_circuit.eval(
            &params,
            &encodings[0],
            &encodings[1..],
            Some(bgg_evaluator),
        );
        // 5. Let `c_f := s^T*A_f + e_{c_f}` in $\mathcal{R}_{q}^{1 \times m}$
        // be the BGG+ encoding corresponding to the output wire of `poly_circuit`.
        let v = ct.c_b_epsilon.concat_rows(&[&result[0].vector]) * fsk.u_f;
        let z = ct.c_u - v.get_row(0)[0].clone();
        z.extract_bits_with_threshold(&params)[0]
    }
}
