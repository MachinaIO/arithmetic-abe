use crate::{
    ciphertext::Ciphertext,
    keys::{FuncSK, MasterPK, MasterSK},
};
use mxx::{
    arithmetic::circuit::ArithmeticCircuit,
    bgg::{
        encoding::BggEncoding,
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    },
    element::PolyElem,
    gadgets::crt::encode_modulo_poly,
    lookup::lwe_eval::LweBggEncodingPltEvaluator,
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
};
use num_bigint::BigUint;
use std::{marker::PhantomData, path::PathBuf, sync::Arc};

const TAG_BGG_PUBKEY: &[u8] = b"BGG_PUBKEY";

pub struct KeyPolicyABE<
    M: PolyMatrix + 'static,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
> {
    pub e_b_sigma: f64,
    pub limb_bit_size: usize,
    pub num_crt_limbs: usize,
    pub crt_depth: usize,
    pub num_eval_slots: usize,
    pub knapsack_size: Option<usize>,
    pub trapdoor_sampler: ST,
    _sh: PhantomData<SH>,
    _su: PhantomData<SU>,
}

impl<
    M: PolyMatrix + 'static,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
> KeyPolicyABE<M, SH, ST, SU>
{
    pub fn new(
        limb_bit_size: usize,
        params: &<M::P as Poly>::Params,
        num_eval_slots: Option<usize>,
        knapsack_size: Option<usize>,
        e_b_sigma: f64,
        trapdoor_sampler: ST,
    ) -> Self {
        let (_, crt_bits, crt_depth) = params.to_crt();
        let num_crt_limbs = crt_bits.div_ceil(limb_bit_size);
        let num_eval_slots = num_eval_slots.unwrap_or(params.ring_dimension() as usize);
        Self {
            limb_bit_size,
            num_crt_limbs,
            crt_depth,
            num_eval_slots,
            knapsack_size,
            e_b_sigma,
            trapdoor_sampler,
            _sh: PhantomData,
            _su: PhantomData,
        }
    }

    pub fn setup(
        &self,
        params: <M::P as Poly>::Params,
        num_inputs: usize,
    ) -> (MasterPK<M>, MasterSK<M, ST>) {
        let seed: [u8; 32] = rand::random();
        let (b_trapdoor, b_matrix) = self.trapdoor_sampler.trapdoor(&params, 1);
        let b_trapdoor = Arc::new(b_trapdoor);
        let b_matrix = Arc::new(b_matrix);
        let uniform_sampler = SU::new();
        let u = uniform_sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let mpk = MasterPK::new(num_inputs, seed, b_matrix, u);
        let msk = MasterSK::new(b_trapdoor);
        (mpk, msk)
    }

    pub fn enc(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        inputs: &[Vec<BigUint>],
        message: &[bool],
    ) -> Ciphertext<M> {
        let num_inputs = inputs.len();
        let uniform_sampler = SU::new();
        let s = uniform_sampler.sample_uniform(&params, 1, 1, DistType::TernaryDist);
        let b_col_size = 2 + params.modulus_digits();
        let c_b_error = {
            let first_part = uniform_sampler.sample_uniform(
                &params,
                1,
                1,
                DistType::GaussDist { sigma: self.e_b_sigma },
            );
            let minus_one = M::P::const_minus_one(&params);
            let second_part = s.clone() * minus_one;
            let third_part = uniform_sampler.sample_uniform(
                &params,
                1,
                b_col_size - 2,
                DistType::GaussDist { sigma: self.e_b_sigma },
            );
            first_part.concat_columns(&[&second_part, &third_part])
        };
        let c_b = s.clone() * mpk.b_matrix.as_ref() + &c_b_error;
        let bgg_encoding_sampler = BGGEncodingSampler::<SU>::new(&params, &s.get_row(0), None);
        let plaintexts = inputs
            .iter()
            .flat_map(|input| {
                assert_eq!(inputs.len(), self.num_eval_slots);
                encode_modulo_poly(self.limb_bit_size, &params, input)
            })
            .collect::<Vec<_>>();
        let num_given_input_polys =
            num_modulo_poly::<M::P>(self.limb_bit_size, &params, num_inputs);
        let reveal_plaintexts = vec![true; num_given_input_polys + 1];
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, 1);
        let pubkeys = bgg_pubkey_sampler.sample(&params, TAG_BGG_PUBKEY, &reveal_plaintexts);
        let bgg_encodings_no_error = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let encode_col_size = params.modulus_digits();
        let knapsack_size = self.knapsack_size.unwrap_or(b_col_size - 1);
        let bgg_encodings = bgg_encodings_no_error
            .into_iter()
            .map(|encode| {
                let mut r_matrix = uniform_sampler.sample_uniform(
                    &params,
                    1,
                    encode_col_size,
                    DistType::TernaryDist,
                );
                r_matrix = r_matrix.concat_rows(&[&M::zero(&params, 1, encode_col_size)]);
                r_matrix = r_matrix.concat_rows(&[&uniform_sampler.sample_uniform(
                    &params,
                    knapsack_size - 1,
                    encode_col_size,
                    DistType::TernaryDist,
                )]);
                if knapsack_size + 1 < b_col_size {
                    r_matrix = r_matrix.concat_rows(&[&M::zero(
                        &params,
                        b_col_size - knapsack_size - 1,
                        encode_col_size,
                    )]);
                }
                let error = c_b_error.clone() * r_matrix;
                let new_vector = encode.vector + error;
                BggEncoding {
                    vector: new_vector,
                    pubkey: encode.pubkey,
                    plaintext: encode.plaintext,
                }
            })
            .collect::<Vec<_>>();
        let ring_dim = params.ring_dimension() as usize;
        assert_eq!(message.len(), self.num_eval_slots, "message length must match num_eval_slots",);
        let mut message_coeffs: Vec<BigUint> =
            message.iter().map(|bit| BigUint::from(*bit as u8)).collect();
        if message_coeffs.len() < ring_dim {
            message_coeffs.resize(ring_dim, BigUint::from(0u8));
        }
        let message_poly = M::P::from_biguints(&params, &message_coeffs);
        let half_q = <M::P as Poly>::Elem::half_q(&params.modulus());
        let half_const = M::P::from_elem_to_constant(&params, &half_q);
        let scaled_message = message_poly * half_const;
        let e_u = uniform_sampler.sample_uniform(
            &params,
            1,
            1,
            DistType::GaussDist { sigma: self.e_b_sigma },
        );
        let c_u = (s.clone() * mpk.u.clone() + e_u).get_row(0)[0].clone() + scaled_message;

        Ciphertext { bgg_encodings, c_b, c_u }
    }

    pub async fn keygen(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, ST>,
        arith_circuit: ArithmeticCircuit<M::P>,
        dir_path: PathBuf,
    ) -> FuncSK<M> {
        let result = arith_circuit
            .evaluate_with_bgg_pubkey::<M, SH, ST, SU>(
                &params,
                mpk.seed,
                dir_path.clone(),
                1,
                mpk.b_matrix.clone(),
                msk.b_trapdoor.clone(),
                self.trapdoor_sampler.clone(),
            )
            .await;
        let a_f = result[0].clone().matrix;
        let u_f = self.trapdoor_sampler.preimage_extend(
            &params,
            &msk.b_trapdoor,
            &mpk.b_matrix,
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
            LweBggEncodingPltEvaluator::<M, SH>::new(mpk.seed, dir_path, ct.c_b.clone());
        let result = arith_circuit.poly_circuit.eval(
            &params,
            &encodings[0],
            &encodings[1..],
            Some(bgg_evaluator),
        );
        // 5. Let `c_f := s^T*A_f + e_{c_f}` in $\mathcal{R}_{q}^{1 \times m}$
        // be the BGG+ encoding corresponding to the output wire of `poly_circuit`.
        let v = ct.c_b.concat_columns(&[&result[0].vector]) * fsk.u_f;
        let z = ct.c_u - &v.get_row(0)[0];
        z.extract_bits_with_threshold(&params)[0]
    }
}

fn num_modulo_poly<P: Poly>(limb_bit_size: usize, params: &P::Params, num_inputs: usize) -> usize {
    let (_, crt_bits, _) = params.to_crt();
    let num_limbs_per_slot = crt_bits.div_ceil(limb_bit_size);
    num_inputs * num_limbs_per_slot
}
