use crate::linear_algebra::{matrix::FMatrix, vector::FVector};

pub struct Orbitals {
    c: FMatrix,  // MO coefficients AO x MO
    e: FVector,  // orbital energies
    nocc: usize, // number of electrons
}
