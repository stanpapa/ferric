use libferric::linear_algebra::{matrix::FMatrix, matrix_container::FMatrixContainer};

pub trait HFSolver {
    fn solve(&mut self, h: &FMatrix, eri: &FMatrixContainer, s: &FMatrix);

    fn guess(&mut self);

    /// Build density as Dμν = \sum_i^{n_occ} Cμi Cνi^T
    fn density(&mut self, s12: &FMatrix);

    fn fock(&mut self, h: &FMatrix, eri: &FMatrixContainer);

    /// Calculate RHF energy as: E0 = 0.5 \sum_{μν} Dμν ( Hμν + Fμν )
    fn energy(&mut self, h: &FMatrix);

    // fn d_rms(&self, d_old: &FMatrix) -> f64;
    fn print_energy(&self, h: &FMatrix);
}
