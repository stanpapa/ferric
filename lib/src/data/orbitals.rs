use crate::linear_algebra::{matrix::FMatrix, vector::FVector};

#[allow(non_camel_case_types)]
pub struct Orbitals {
    C: [FMatrix; 2],  // MO coefficients AO x MO
    E: [FVector; 2],  // orbital energies
    ON: [FVector; 2], // occupation number vectors

    // number of electrons
    N: usize,         // total number of electrons
    Nα: usize,       // number of alpha electrons
    Nβ: usize,       // number of electrons
    homo: [usize; 2], // number of electrons

    num_op: usize,
}

// impl Orbitals {
//     pub fn new(guess: &FMatrix) -> Self {}
// }
