use crate::{
    ciphertext::Ciphertext,
    keys::{FuncSK, MasterPK, MasterSK},
};
use log::info;
use mxx::{
    // arithmetic::circuit::ArithmeticCircuit,
    bgg::{
        encoding::BggEncoding,
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    },
    circuit::PolyCircuit,
    element::PolyElem,
    gadgets::arith::nested_rns::{NestedRnsPoly, NestedRnsPolyContext, encode_nested_rns_poly},
    lookup::ggh15_eval::{GGH15BGGEncodingPltEvaluator, GGH15BGGPubKeyPltEvaluator},
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
    storage::write::{init_storage_system, wait_for_all_writes},
};
use num_bigint::BigUint;
use std::{marker::PhantomData, path::PathBuf, sync::Arc};

const TAG_BGG_PUBKEY: &[u8] = b"BGG_PUBKEY";

pub struct KeyPolicyABE<
    M: PolyMatrix + 'static,
    HS: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    TS: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    US: PolyUniformSampler<M = M> + Send + Sync,
> {
    pub p_moduli_bits: usize,
    // pub p_moduli_depth: usize,
    pub scale: u64,
    pub crt_depth: usize,
    pub knapsack_size: Option<usize>,
    pub trapdoor_sigma: f64,
    pub e_b_sigma: f64,
    // pub trapdoor_sampler: TS,
    _hs: PhantomData<HS>,
    _us: PhantomData<US>,
    _ts: PhantomData<TS>,
}

impl<
    M: PolyMatrix + 'static,
    HS: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    TS: PolyTrapdoorSampler<M = M> + Clone + Send + Sync,
    US: PolyUniformSampler<M = M> + Send + Sync,
> KeyPolicyABE<M, HS, TS, US>
{
    pub fn new(
        p_moduli_bits: usize,
        scale: u64,
        params: &<M::P as Poly>::Params,
        knapsack_size: Option<usize>,
        trapdoor_sigma: f64,
        e_b_sigma: f64,
    ) -> Self {
        assert!(p_moduli_bits > 1, "p_moduli_bits must be at least 2 for NestedRns");
        let (_, _, crt_depth) = params.to_crt();
        // NestedRns chooses the number of level-1 moduli as ceil(2 * |q_bits| / (p_bits - 1))
        // let p_moduli_depth = (2 * crt_bits).div_ceil(p_moduli_bits - 1);
        Self {
            p_moduli_bits,
            // p_moduli_depth,
            scale,
            crt_depth,
            knapsack_size,
            trapdoor_sigma,
            e_b_sigma,
            _hs: PhantomData,
            _us: PhantomData,
            _ts: PhantomData,
        }
    }

    pub fn setup(
        &self,
        params: &<M::P as Poly>::Params,
        num_inputs: usize,
    ) -> (MasterPK<M>, MasterSK<M, TS>) {
        let seed: [u8; 32] = rand::random();
        let trapdoor_sampler = TS::new(params, self.trapdoor_sigma);
        let (b_trapdoor, b_matrix) = trapdoor_sampler.trapdoor(params, 1);
        let b_trapdoor = Arc::new(b_trapdoor);
        let b_matrix = Arc::new(b_matrix);
        let uniform_sampler = US::new();
        let u = uniform_sampler.sample_uniform(params, 1, 1, DistType::FinRingDist);
        let mpk = MasterPK::new(num_inputs, seed, b_matrix, u);
        let msk = MasterSK::new(b_trapdoor);
        (mpk, msk)
    }

    pub fn enc(
        &self,
        params: &<M::P as Poly>::Params,
        mpk: MasterPK<M>,
        inputs: &[BigUint],
        message: bool,
    ) -> Ciphertext<M> {
        let num_inputs = inputs.len();
        assert_eq!(
            num_inputs, mpk.num_inputs,
            "provided inputs ({num_inputs}) must match mpk.num_inputs ({})",
            mpk.num_inputs
        );
        let uniform_sampler = US::new();
        let s = uniform_sampler.sample_uniform(&params, 1, 1, DistType::TernaryDist);
        let b_col_size = 2 + params.modulus_digits();
        let c_b_error: M = {
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
        let bgg_encoding_sampler = BGGEncodingSampler::<US>::new(params, &s.get_row(0), None);
        // let (_, _, crt_depth) = params.to_crt();
        // let p_moduli_depth = (2 * crt_bits).div_ceil(self.p_moduli_bits - 1);
        let plaintexts = inputs
            .iter()
            .flat_map(|input| encode_nested_rns_poly(self.p_moduli_bits, params, input))
            .collect::<Vec<_>>();
        // let expected_plaintexts = mpk.num_inputs * crt_depth * self.p_moduli_depth;
        // assert_eq!(
        //     plaintexts.len(),
        //     expected_plaintexts,
        //     "plaintext count ({}) must equal num_inputs * crt_depth * p_moduli_depth ({})",
        //     plaintexts.len(),
        //     expected_plaintexts
        // );
        let reveal_plaintexts = vec![true; plaintexts.len()];
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, HS>::new(mpk.seed, 1);
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
        // let ring_dim = params.ring_dimension() as usize;
        // let mut message_coeffs: Vec<BigUint> =
        //     message.iter().map(|bit| BigUint::from(*bit as u8)).collect();
        // if message_coeffs.len() < ring_dim {
        //     message_coeffs.resize(ring_dim, BigUint::from(0u8));
        // }
        let message_poly = M::P::from_usize_to_constant(&params, message as usize);
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
        params: &<M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, TS>,
        height: u32,
        dir_path: PathBuf,
    ) -> FuncSK<M> {
        init_storage_system();
        let circuit = {
            let mut circuit = PolyCircuit::<M::P>::new();
            let ctx = Arc::new(NestedRnsPolyContext::setup(
                &mut circuit,
                &params,
                self.p_moduli_bits,
                self.scale,
                false,
            ));
            NestedRnsPoly::benchmark_multiplication_tree(
                ctx,
                params,
                &mut circuit,
                height as usize,
            );
            circuit
        };
        let plt_evaluator = GGH15BGGPubKeyPltEvaluator::<M, US, HS, TS>::new(
            mpk.seed,
            self.trapdoor_sigma,
            self.e_b_sigma,
            &params,
            mpk.b_matrix.clone(),
            msk.b_trapdoor.clone(),
            dir_path.clone(),
            false,
        );
        let reveal_plaintexts = vec![true; circuit.num_input()];
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, HS>::new(mpk.seed, 1);
        let pubkeys = bgg_pubkey_sampler.sample(&params, TAG_BGG_PUBKEY, &reveal_plaintexts);
        let result = circuit.eval(params, &pubkeys[0], &pubkeys[1..], Some(plt_evaluator));
        info!("finished evaluation of pubkeys");
        wait_for_all_writes(dir_path.clone()).await.unwrap();
        info!("finished write files");

        let a_f = result[0].clone().matrix;
        let trapdoor_sampler = TS::new(params, self.trapdoor_sigma);
        let u_f =
            trapdoor_sampler.preimage_extend(&params, &msk.b_trapdoor, &mpk.b_matrix, &a_f, &mpk.u);
        assert_eq!(result.len(), 1);
        FuncSK { a_f, u_f, dir_path }
    }

    pub fn dec(
        &self,
        params: &<M::P as Poly>::Params,
        ct: Ciphertext<M>,
        mpk: MasterPK<M>,
        fsk: FuncSK<M>,
        height: u32,
    ) -> bool {
        init_storage_system();
        let circuit = {
            let mut circuit = PolyCircuit::<M::P>::new();
            let ctx = Arc::new(NestedRnsPolyContext::setup(
                &mut circuit,
                &params,
                self.p_moduli_bits,
                self.scale,
                false,
            ));
            NestedRnsPoly::benchmark_multiplication_tree(
                ctx,
                params,
                &mut circuit,
                height as usize,
            );
            circuit
        };
        let encodings = &ct.bgg_encodings[..];
        assert_eq!(
            encodings.len(),
            circuit.num_input() + 1,
            "ciphertext must contain exactly 1 + circuit.num_input() encodings"
        );
        let dir_path: PathBuf = fsk.dir_path;
        let bgg_evaluator = GGH15BGGEncodingPltEvaluator::<M, HS>::new(
            mpk.seed,
            &params,
            dir_path,
            mpk.b_matrix.row_size(),
            ct.c_b.clone(),
        );
        let result = circuit.eval(params, &encodings[0], &encodings[1..], Some(bgg_evaluator));
        // 5. Let `c_f := s^T*A_f + e_{c_f}` in $\mathcal{R}_{q}^{1 \times m}$
        // be the BGG+ encoding corresponding to the output wire of `poly_circuit`.
        let v = ct.c_b.concat_columns(&[&result[0].vector]) * fsk.u_f;
        let z = ct.c_u - &v.get_row(0)[0];
        z.extract_bits_with_threshold(params)[0]
    }
}
