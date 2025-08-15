use crate::{
    ciphertext::Ciphertext,
    circuit::ArithmeticCircuit,
    functional_key::FuncSK,
    master_key::{MasterPK, MasterSK},
};
use mxx::circuit::gate::GateId;
use mxx::element::PolyElem;
use mxx::lookup::simple_eval::SimpleBggPubKeyEvaluator;
use mxx::utils::log_mem;
use mxx::{
    bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    circuit::PolyCircuit,
    gadgets::{
        crt::{CrtContext, CrtPoly},
        lt_isolate::LtIsolateGadget,
    },
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
};
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
        let (b_epsilon_trapdoor, b_epsilon) = self.trapdoor_sampler.trapdoor(&params, self.d);
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
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs.div_ceil(mpk.packed_limbs);
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let s = &self
            .uniform_sampler
            .sample_uniform(&params, 1, self.d, DistType::BitDist);
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);

        // ===== CRT =====
        let mut outputs: Vec<GateId> = Vec::with_capacity(mpk.num_inputs);
        let mut circuit = PolyCircuit::<M::P>::new();
        let inputs_gates = circuit.input(total_limbs);
        let ctx = Arc::new(CrtContext::setup(&mut circuit, &params, self.limb_bit_size));
        let num_crt_limbs = self.num_crt_limbs;
        for i in 0..mpk.num_inputs {
            let crt_poly = CrtPoly::from_inputs_interleaved(
                &mut circuit,
                ctx.clone(),
                &inputs_gates,
                num_crt_limbs,
                i,
                mpk.num_inputs,
            );
            outputs.extend(crt_poly.limb());
        }
        assert_eq!(outputs.len(), total_limbs);
        circuit.output(outputs);
        let packed_inputs: Vec<M::P> = CrtPoly::<M::P>::generate_input_values_from_single(
            &params,
            &ctx,
            inputs,
            self.limb_bit_size,
        );
        assert_eq!(num_packed_poly_inputs, packed_inputs.len());

        // ===== Enc =====
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
        let bgg_sampler = BGGEncodingSampler::new(&params, &s.get_row(0), SU::new(), self.p_sigma);
        let bgg_encodings = bgg_sampler.sample(&params, &pubkeys, &packed_inputs);
        let c_b_epsilon = s.clone() * mpk.b_epsilon + c_b_epsilon_error;
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

    pub fn keygen(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, ST>,
        mut arith_circuit: ArithmeticCircuit<M::P>,
    ) -> FuncSK<M> {
        let ring_dim = params.ring_dimension() as usize;
        let k = arith_circuit.packed_limbs.saturating_sub(1);
        let lt_isolate_id = LtIsolateGadget::register_general_lt_isolate_lookup(
            &mut arith_circuit.original_circuit,
            &params,
            k,
        );
        arith_circuit.to_poly_circuit(lt_isolate_id, ring_dim);
        log_mem("finish to_poly_circuit");
        let poly_circuit = arith_circuit.original_circuit.clone();
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs.div_ceil(mpk.packed_limbs);
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        let b_epsilon_trapdoor = Arc::new(msk.b_epsilon_trapdoor);
        let b_epsilon = Arc::new(mpk.b_epsilon);
        let dir_path: PathBuf = "keygen".into();
        let bgg_plt_evaluator = SimpleBggPubKeyEvaluator::<M, SH, SU, ST>::new(
            mpk.seed,
            self.trapdoor_sampler.clone(),
            b_epsilon.clone(),
            b_epsilon_trapdoor.clone(),
            dir_path.clone(),
        );
        // todo: error on
        // assertion `left == right` failed: number of provided inputs must match circuit inputs left: 3 right: 60
        // basically sampled pubkeys is treating as integer arithmetic circuit input also packed, poly_circuit is still accepting as crt + limb decomposed form.
        let result =
            poly_circuit.eval(&params, &pubkeys[0], &pubkeys[1..], Some(bgg_plt_evaluator));
        assert_eq!(result.len(), 1);
        let a_f = result[0].clone().matrix;
        let u_f = self.trapdoor_sampler.preimage_extend(
            &params,
            &b_epsilon_trapdoor,
            &b_epsilon,
            &a_f,
            &mpk.u,
        );

        FuncSK { a_f, u_f, dir_path }
    }

    pub fn dec(
        &self,
        params: <M::P as Poly>::Params,
        ct: Ciphertext<M>,
        mpk: MasterPK<M>,
        fsk: FuncSK<M>,
        mut arith_circuit: ArithmeticCircuit<M::P>,
    ) -> bool {
        let ring_dim = params.ring_dimension() as usize;
        let k = arith_circuit.packed_limbs.saturating_sub(1);
        let lt_isolate_id = LtIsolateGadget::register_general_lt_isolate_lookup(
            &mut arith_circuit.original_circuit,
            &params,
            k,
        );
        arith_circuit.to_poly_circuit(lt_isolate_id, ring_dim);
        let poly_circuit = arith_circuit.original_circuit.clone();
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs.div_ceil(mpk.packed_limbs);
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        // TODO: General question on this step about from this, step 3 we get "bgg public key", and so it refer step 4 is evaluating over public keys. and then suddenly step 5 is refering I can get BGG+ encoding from output wire of circuit which is conflicting from my understandation. Did i missed smth?
        // let bgg_plt_evaluator = SimpleBggPubKeyEvaluator::<M, SH, SU, ST>::new(
        //     mpk.seed,
        //     self.trapdoor_sampler,
        //     Arc::new(mpk.b_epsilon),
        //     Arc::new(msk.b_epsilon_trapdoor),
        //     "keygen".into(),
        // );
        let result = poly_circuit.eval(
            &params,
            &pubkeys[0],
            &pubkeys[1..],
            None::<SimpleBggPubKeyEvaluator<M, SH, SU, ST>>,
        );
        // 5. Let `c_f := s^T*A_f + e_{c_f}` in $\mathcal{R}_{q}^{1 \times m}$
        // be the BGG+ encoding corresponding to the output wire of `poly_circuit`.
        // TODO: Compute `v := (c_{B_{\epsilon}}, c_f) * u_f = s^{T} * u + e_{c_f} * u_f`. <--- how do i concat Also how do i get c_f esp i'm evaluating over BggPubkey
        let v = ct.c_b_epsilon.concat_rows(&[&result[0].matrix]) * fsk.u_f;
        let z = ct.c_u - v.get_row(0)[0].clone();
        z.extract_bits_with_threshold(&params)[0]
    }
}
