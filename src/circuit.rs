use mxx::{circuit::PolyCircuit, gadgets::lt_isolate::LtIsolateGadget, poly::Poly};

#[derive(Clone)]
pub struct ArithmeticCircuit<P: Poly> {
    pub limb_bit_size: usize,
    pub crt_bits: usize,
    pub crt_depth: usize,
    pub num_crt_limbs: usize,
    pub packed_limbs: usize,
    pub num_inputs: usize,
    pub original_circuit: PolyCircuit<P>,
}

impl<P: Poly> ArithmeticCircuit<P> {
    pub fn to_poly_circuit(&mut self, lt_isolate_id: usize, ring_dim: usize) {
        let k = self.packed_limbs.saturating_sub(1);
        let input_gates = self.original_circuit.input(self.num_inputs);

        let mut decomposed_outputs = Vec::with_capacity(input_gates.len() * self.packed_limbs);
        for g in input_gates {
            let isolated = LtIsolateGadget::isolate_coeffs(
                &mut self.original_circuit,
                g,
                k,
                lt_isolate_id,
                ring_dim,
            );
            decomposed_outputs.extend(isolated);
        }
        self.original_circuit.output(decomposed_outputs);
    }

    pub fn create_isolate_coeffs_subcircuit(
        &self,
        params: &<P as Poly>::Params,
        k: usize,
        ring_dim: usize,
    ) -> PolyCircuit<P> {
        let mut sc = PolyCircuit::<P>::new();
        let lt_isolate_id = LtIsolateGadget::register_general_lt_isolate_lookup(&mut sc, params, k);
        let a_poly = sc.input(1)[0];
        let isolated = LtIsolateGadget::isolate_coeffs(&mut sc, a_poly, k, lt_isolate_id, ring_dim);
        sc.output(isolated);
        sc
    }
}
