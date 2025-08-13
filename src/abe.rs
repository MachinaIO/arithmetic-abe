use mxx::element::PolyElem;
use mxx::{
    bgg::{
        encoding::BggEncoding,
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    },
    circuit::PolyCircuit,
    gadgets::crt::{CrtContext, CrtPoly},
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
    utils::create_bit_random_poly,
};
use std::{marker::PhantomData, sync::Arc};

use crate::{
    ciphertext::Ciphertext,
    circuit::ArithmeticCircuit,
    functional_key::FuncSK,
    master_key::{MasterPK, MasterSK},
};

const TAG_BGG_PUBKEY: &[u8] = b"BGG_PUBKEY";

pub struct KeyPolicyABE<
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M>,
    ST: PolyTrapdoorSampler<M = M>,
    SU: PolyUniformSampler<M = M>,
> {
    limb_bit_size: usize,
    num_crt_limbs: usize,
    crt_depth: usize,
    packed_limb: usize,
    d: usize,
    hash_sampler: SH,
    trapdoor_sampler: ST,
    uniform_sampler: SU,
    _i: PhantomData<M>,
}

impl<
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M>,
    ST: PolyTrapdoorSampler<M = M>,
    SU: PolyUniformSampler<M = M> + Copy,
> KeyPolicyABE<M, SH, ST, SU>
{
    fn setup(
        &self,
        params: <M::P as Poly>::Params,
        num_inputs: usize,
        packed_limbs: usize,
    ) -> (MasterPK<M>, MasterSK<M, ST>) {
        let seed: [u8; 32] = rand::random();
        let num_packed_poly_inputs =
            (num_inputs * self.num_crt_limbs * self.crt_depth) / packed_limbs;
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(seed, self.d);
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        let (b_epsilon_trapdoor, b_epsilon) = self.trapdoor_sampler.trapdoor(&params, self.d);
        let u = self
            .uniform_sampler
            .sample_uniform(&params, self.d, 1, DistType::BitDist);
        let mpk = MasterPK::new(num_inputs, packed_limbs, seed, b_epsilon, u);
        let msk = MasterSK::new(b_epsilon_trapdoor);
        (mpk, msk)
    }

    fn enc(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        inputs: &[<M::P as Poly>::Elem],
        message: bool,
        p_sigma: f64,
    ) -> Ciphertext<M> {
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs / mpk.packed_limbs;
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let s = self
            .uniform_sampler
            .sample_uniform(&params, self.d, 1, DistType::BitDist);
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        let (moduli, crt_bits, _) = params.to_crt();
        debug_assert_eq!(moduli.len(), self.crt_depth);
        // 4. For every `i in 0..num_packed_poly_inputs`, decompose `inputs[i]`
        // into the CRT representation + the big integer representation,
        // stored in `crt_inputs[i*(num_crt_limbs*crt_depth)..(i+1)*(num_crt_limbs*crt_depth)]`.
        let mut crt_inputs: Vec<CrtPoly<M::P>> = Vec::with_capacity(num_packed_poly_inputs);
        let mut circuit = PolyCircuit::<M::P>::new();
        let inputs = circuit.input(total_limbs);
        let ctx = Arc::new(CrtContext::setup(&mut circuit, &params, self.limb_bit_size));
        for i in 0..num_packed_poly_inputs {
            // let target_input = inputs[i].value();
            let crt_poly = CrtPoly::from_inputs_interleaved(
                &mut circuit,
                ctx.clone(),
                &inputs,
                self.num_crt_limbs,
                i,
                num_packed_poly_inputs,
            );
            crt_inputs.push(crt_poly);
        }
        // todo: 5. For every `i in 0..num_packed_poly_inputs`, make a packed polynomial `packed_inputs[i]` from the `packed_limbs` integers in `crt_inputs`.
        let mut packed_inputs: Vec<M::P> = Vec::with_capacity(num_packed_poly_inputs);

        let s_bar = self
            .uniform_sampler
            .sample_uniform(&params, 1, self.d, DistType::BitDist);
        // todo: gauss_sigma and p_sigma
        let bgg_sampler =
            BGGEncodingSampler::new(&params, &s_bar.get_row(0), self.uniform_sampler, p_sigma);
        let bgg_encodings = bgg_sampler.sample(&params, &pubkeys, &packed_inputs);
        let c_b_epsilon_error = self.uniform_sampler.sample_uniform(
            &params,
            1,
            self.d * (2 + params.modulus_digits()),
            DistType::GaussDist { sigma: p_sigma },
        );
        let c_b_epsilon = s_bar.clone() * mpk.b_epsilon + c_b_epsilon_error;
        let boolean_msg = if message {
            <M::P as Poly>::Elem::one(&params.modulus())
        } else {
            <M::P as Poly>::Elem::zero(&params.modulus())
        };
        let scale = M::P::from_elem_to_constant(
            &params,
            &(<M::P as Poly>::Elem::half_q(&params.modulus()) * boolean_msg),
        );
        let c_u = (s_bar * mpk.u).get_row(0)[0].clone() + scale;

        Ciphertext {
            bgg_encodings,
            c_b_epsilon,
            c_u,
        }
    }

    fn keygen(
        &self,
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, ST>,
        circuit: ArithmeticCircuit,
    ) -> FuncSK {
        // 1. Convert `arith_circuit` into `poly_circuit: PolyCircuit` in the way described above.
        // 2. Reconstruct `A_{1}, A_{x_{0}}, \dots, A_{x_{num_packed_poly_inputs-1}}` from `seed` with the hash sampler.
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs / mpk.packed_limbs;
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        // 3. Evaluate `poly_circuit` on the above BGG+ public matrices. During this process, required preimages and matrices for every lookup gate should be generated.
        // 4. Let `A_f` in $\mathcal{R}_{q}^{d \times m}$ be the BGG+ public matrix
        // corresponding to the output wire of `poly_circuit`.
        // Generate a preimage `u_f := [B_{\epsilon},A_f]^{-1}(u)`
        // by calling the [`preimage_extend` function](https://github.com/MachinaIO/mxx/blob/main/src/sampler/trapdoor/sampler.rs#L159C8-L159C23).
        // 5. Output the preimages/matrices generated in Step 3 along with `arith_circuit` and `A_f` as `fsk: FuncSk`.

        todo!()
    }

    fn dec(&self, params: <M::P as Poly>::Params, mpk: MasterPK<M>, fsk: FuncSK) -> bool {
        // 2. Convert `arith_circuit` into `poly_circuit: PolyCircuit` in the way described above.
        // 3. Reconstruct `A_{1}, A_{x_{0}}, \dots, A_{x_{num_packed_poly_inputs-1}}` from `seed` with the hash sampler.
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let total_limbs = self.num_crt_limbs * self.crt_depth * mpk.num_inputs;
        let num_packed_poly_inputs = total_limbs / mpk.packed_limbs;
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        // 4. Evaluate `poly_circuit` on the above BGG+ public matrices.
        // 5. Let `c_f := s^T*A_f + e_{c_f}` in $\mathcal{R}_{q}^{1 \times m}$
        // be the BGG+ encoding corresponding to the output wire of `poly_circuit`.
        // Compute `v := (c_{B_{\epsilon}}, c_f) * u_f = s^{T} * u + e_{c_f} * u_f`.
        // 6. Compute `z = c_{u} - v = e_{c_{u}} - e_{c_f} * u_f + \frac{q}{2}\mu`.
        // 7. Output `\mu=1` if the constant term of `z` is larger than the threshold and `\mu=0` otherwise.

        todo!()
    }
}
