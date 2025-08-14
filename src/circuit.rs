use mxx::{circuit::PolyCircuit, poly::Poly};

pub struct ArithmeticCircuit<P: Poly> {
    limb_bit_size: usize,
    crt_bits: usize,
    crt_depth: usize,
    num_crt_limbs: usize,
    pub packed_limbs: usize,
    num_inputs: usize,
    pub original_circuit: PolyCircuit<P>,
}

impl<P: Poly> ArithmeticCircuit<P> {
    pub fn to_poly_circuit(mut self, lt_isolate_id: usize, ring_dim: usize) -> PolyCircuit<P> {
        let k = self.packed_limbs.saturating_sub(1);
        let input_gates = self.original_circuit.input(self.num_inputs);

        let mut decomposed_outputs = Vec::with_capacity(input_gates.len() * self.packed_limbs);
        for g in input_gates {
            let isolated = self
                .original_circuit
                .isolate_coeffs(g, k, lt_isolate_id, ring_dim);
            decomposed_outputs.extend(isolated);
        }
        self.original_circuit.output(decomposed_outputs);
        self.original_circuit
    }

    // If you *need* a subcircuit, re-register the lookup *inside* it.
    pub fn create_isolate_coeffs_subcircuit(
        &self,
        params: &<P as Poly>::Params,
        k: usize,        // usually self.packed_limbs - 1
        ring_dim: usize, // params.ring_dimension() as usize
    ) -> PolyCircuit<P> {
        let mut sc = PolyCircuit::<P>::new();
        let lt_isolate_id = sc.register_general_lt_isolate_lookup(params, k);
        let a_poly = sc.input(1)[0];
        let isolated = sc.isolate_coeffs(a_poly, k, lt_isolate_id, ring_dim);
        sc.output(isolated);
        sc
    }
}
