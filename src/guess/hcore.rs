use libferric::{
    geometry::Geometry,
    gto_integrals::one_electron::OneElectronKernel,
    linear_algebra::{diagonalize::DiagonalizeSym, matrix::FMatrix, power::PowerSym},
    HFType,
};

pub fn guess(basename: &str, hf: &HFType, geometry: &Geometry) {
    print!("Loading integrals  ... ");
    let h = FMatrix::retrieve(OneElectronKernel::HCore.to_filename());
    let s = FMatrix::retrieve(OneElectronKernel::Overlap.to_filename());
    let s12 = s.powf_sym(-0.5);
    println!("done");

    print!("Construction guess ... ");
    match hf {
        HFType::RHF => guess_rhf(basename, &h, &s12, &geometry.n_electrons),
        HFType::UHF => guess_uhf(
            basename,
            &h,
            &s12,
            &[geometry.n_electrons_alpha, geometry.n_electrons_beta],
        ),
        _ => panic!("Guess: unsupported HFType for guess"),
    }
    print!("done");
}

fn guess_rhf(basename: &str, h: &FMatrix, s12: &FMatrix, nel: &usize) {
    let homo = nel / 2;

    // --------------------------------
    // build initial guess density
    // --------------------------------
    // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
    // S-1/2 = symmetric
    let f = s12 * (h * s12);

    // diagonalize F': C'^T F' C' = eps
    let (_, cprime) = f.diagonalize_sym();
    let c = s12 * cprime;

    // --------------------------------
    // Initial density
    // --------------------------------
    let c_occ = c.slice(0, c.rows - 1, 0, homo - 1);
    let d = &c_occ * c_occ.transposed();

    // --------------------------------
    // Store guess on disk
    // --------------------------------
    d.store(&format!("{basename}.p.tmp"));
}

fn guess_uhf(basename: &str, h: &FMatrix, s12: &FMatrix, homo: &[usize; 2]) {
    (0..2).for_each(|op| {
        // --------------------------------
        // build initial guess density
        // --------------------------------
        // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
        // S-1/2 = symmetric
        let f = s12 * (h * s12);

        // diagonalize F': C'^T F' C' = eps
        let (_, cprime) = f.diagonalize_sym();
        let c = s12 * cprime;

        // --------------------------------
        // Initial density
        // --------------------------------
        let c_occ = c.slice(0, c.rows - 1, 0, homo[op] - 1);
        let d = &c_occ * c_occ.transposed();

        // --------------------------------
        // Store guess on disk
        // --------------------------------
        d.store(&format!("{basename}.p{op}.tmp",));
    });
}
