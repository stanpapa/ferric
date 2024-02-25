mod diis;
mod fock;
mod rhf;
mod solver;
mod uhf;

use crate::rhf::RHFSolver;
use crate::solver::HFSolver;
use crate::uhf::UHFSolver;

use libferric::{
    geometry::{molecule::Molecule, Geometry},
    gto_basis_sets::{basis::Basis, load_basis_set, BasisSet},
    gto_integrals::{
        integral_interface::IntegralInterface, one_electron::OneElectronKernel,
        two_electron::TwoElectronKernel,
    },
    linear_algebra::{
        matrix::FMatrix, matrix_container::FMatrixContainer, matrix_container::FMatrixSymContainer,
        matrix_symmetric::FMatrixSym, power::Power,
    },
    HFType,
    HFType::{RHF, UHF},
};

use std::error;

fn main() -> Result<(), Box<dyn error::Error>> {
    println!("SCF Module");
    let basename = "input";

    // hard-coded base name for now
    let geometry = Geometry::retrieve(basename);

    let basis = Basis::retrieve(basename); // todo!()
    let integrals = IntegralInterface::new(&basis, geometry.molecule.atoms());

    // read integrals from disk
    let h_core = FMatrixSym::retrieve(OneElectronKernel::HCore.to_filename());
    let h = FMatrix::from(h_core); // FMatrixSym -> FMatrix
    let s = FMatrixSym::retrieve(OneElectronKernel::Overlap.to_filename());
    let eri = integrals.calc_two_electron_integral(TwoElectronKernel::ERI); // todo: verify

    let hf = RHF;
    println!("{} calculation", hf);

    match hf {
        RHF => {
            let mut rhf = RHFSolver::new(&h, &geometry);
            rhf.solve(&h, &eri, &s);
        }
        UHF => {
            let mut uhf = UHFSolver::new(&[h.clone(), h.clone()], &geometry);
            uhf.solve(&h, &eri, &s);
        }
        _ => panic!("{} not implemented", hf),
    }

    Ok(())
}
