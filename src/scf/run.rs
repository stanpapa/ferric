use super::{input::SCFInput, rhf::RHFSolver, solver::HFSolver, uhf::UHFSolver};

use libferric::{
    geometry::Geometry,
    gto_basis_sets::basis::Basis,
    gto_integrals::{
        integral_interface::IntegralInterface, one_electron::OneElectronKernel,
        two_electron::TwoElectronKernel,
    },
    linear_algebra::{matrix::FMatrix, matrix_symmetric::FMatrixSym},
    HFType::{RHF, UHF},
};

use std::error;

pub fn run(basename: &str, scf_input: SCFInput) -> Result<(), Box<dyn error::Error>> {
    println!("SCF Module");

    // read data
    let geometry = Geometry::retrieve(basename);
    let basis = Basis::retrieve(basename);
    let integrals = IntegralInterface::new(&basis, geometry.molecule.atoms());

    // read integrals from disk
    let h_core = FMatrixSym::retrieve(OneElectronKernel::HCore.to_filename());
    let h = FMatrix::from(h_core); // FMatrixSym -> FMatrix
    let s = FMatrixSym::retrieve(OneElectronKernel::Overlap.to_filename());
    let eri = integrals.calc_two_electron_integral(TwoElectronKernel::ERI); // todo: verify

    println!("{} calculation", scf_input.hf);

    let mut solver = set_solver(scf_input.clone(), &h, &geometry);
    solver.solve(&h, &eri, &s);

    Ok(())
}

fn set_solver(scf_input: SCFInput, h: &FMatrix, geometry: &Geometry) -> Box<dyn HFSolver> {
    match scf_input.hf {
        RHF => Box::new(RHFSolver::new(h, geometry, scf_input)),
        UHF => Box::new(UHFSolver::new(&[h.clone(), h.clone()], geometry, scf_input)),
        _ => panic!("{} not implemented", scf_input.hf),
    }
}
