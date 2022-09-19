use crate::{
    geometry::atom::Atom,
    gto_integrals::{
        kinetic_energy::kinetic_energy, nuclear_electron_attraction::nuclear_electron_attraction,
    },
};

pub fn h_core(
    a: &f64,
    ml_a: &[i16; 3],
    a_origin: &[f64; 3],
    b: &f64,
    ml_b: &[i16; 3],
    b_origin: &[f64; 3],
    atoms: &[Atom],
) -> f64 {
    let kinetic = kinetic_energy(a, ml_a, a_origin, b, ml_b, b_origin);
    let nuclear_attraction =
        nuclear_electron_attraction(a, ml_a, a_origin, b, ml_b, b_origin, atoms);

    kinetic + nuclear_attraction
}
