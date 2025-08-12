use mxx::{
    bgg::sampler::BGGPublicKeySampler,
    matrix::PolyMatrix,
    poly::{Poly, PolyParams},
    sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
};
use std::marker::PhantomData;

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
    SU: PolyUniformSampler<M = M>,
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
    ) -> Ciphertext {
        let num_packed_poly_inputs =
            (mpk.num_inputs * self.num_crt_limbs * self.crt_depth) / mpk.packed_limbs;
        let reveal_plaintexts = vec![true; num_packed_poly_inputs + 1];
        let s = self
            .uniform_sampler
            .sample_uniform(&params, self.d, 1, DistType::BitDist);
        let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(mpk.seed, self.d);
        let pubkeys = bgg_pubkey_sampler.sample(&params, &TAG_BGG_PUBKEY, &reveal_plaintexts);
        // 4. For every `i in 0..num_packed_poly_inputs`, decompose `inputs[i]` into the CRT representation + the big integer representation, stored in `crt_inputs[i*(num_crt_limbs*crt_depth)..(i+1)*(num_crt_limbs*crt_depth)]`.
        // 5. For every `i in 0..num_packed_poly_inputs`, make a packed polynomial `packed_inputs[i]` from the `packed_limbs` integers in `crt_inputs`.
        // 6. Construct BGG+ encodings: `c_{1}:=s^T*(A_{1}-gadget) + e_{c_1}, c_{x_1}:=s^T * (A_{x_1} - packed_inputs[1]*gadget) + e_{c_{x_1}}, ..., c_{x_{num_packed_poly_inputs-1}}:=s^T * (A_{x_{num_packed_poly_inputs-1}} - packed_inputs[1]*gadget) + e_{c_{x_{num_packed_poly_inputs-1}}}`, which can be sampled through `BGGEncodingSampler`.
        // 7. Compute `c_{b_{\epsilon}} := s^{T} * B_{\epsilon} + e_{c_{b_{\epsilon}}}`.
        // 8. Compute `c_{u} := s^{T} * u + e_{c_{u}} + \frac{q}{2}\mu`.
        // 9. Output `ct: Ciphertext := {c_{1}, c_{x_1}, \dots, c_{x_{num_packed_poly_inputs-1}}, c_{b_{\epsilon}}, c_{u}}`.

        todo!()
    }

    fn keygen(
        params: <M::P as Poly>::Params,
        mpk: MasterPK<M>,
        msk: MasterSK<M, ST>,
        circuit: ArithmeticCircuit,
    ) -> FuncSK {
        todo!()
    }

    fn dec(params: <M::P as Poly>::Params, mpk: MasterPK<M>, fsk: FuncSK) -> bool {
        todo!()
    }
}
