use ferric_lib::{
    geometry::molecule::Molecule,
    linear_algebra::{
        matrix::FMatrix, matrix_container::FMatrixSymContainer, matrix_symmetric::FMatrixSym,
    },
};

pub trait HFSolver {
    fn solve(&mut self, h: &FMatrix, eri: &FMatrixSymContainer, s: &FMatrixSym);

    fn guess(&mut self, h: &FMatrix, s12: &FMatrix);

    /// Build density as Dμν = \sum_i^{n_occ} Cμi Cνi^T
    fn density(&mut self);

    fn fock(&mut self, h: &FMatrix, eri: &FMatrixSymContainer);

    /// Calculate RHF energy as: E0 = 0.5 \sum_{μν} Dμν ( Hμν + Fμν )
    fn energy(&mut self, h: &FMatrix);

    // fn d_rms(&self, d_old: &FMatrix) -> f64;
    fn print_energy(&self, h: &FMatrix);
}
