use super::{input::SCFInput, rhf::RHFSolver, solver::HFSolver, uhf::UHFSolver};

use libferric::{
    geometry::Geometry,
    gto_basis_sets::basis::Basis,
    gto_integrals::{
        integral_interface::IntegralInterface, one_electron::OneElectronKernel,
        two_electron::TwoElectronKernel,
    },
    linear_algebra::matrix::FMatrix,
    HFType::{RHF, UHF},
};

use std::{error, time::Instant};

pub fn run(basename: &str, scf_input: SCFInput) -> Result<(), Box<dyn error::Error>> {
    println!("=====================================");
    println!("               SCF");
    println!("=====================================\n");

    // read data
    let geometry = Geometry::retrieve(basename);
    let basis = Basis::retrieve(basename);
    let integrals = IntegralInterface::new(&basis, geometry.molecule.atoms());

    println!("-------------------------------------");
    println!("             Settings");
    println!("-------------------------------------\n");
    println!("General\n-------");
    println!("Hartree-Fock Type:                {}", scf_input.hf);
    println!("Total Charge:                     {}", geometry.charge);
    println!(
        "Multiplicity:                     {}",
        geometry.multiplicity
    );
    match scf_input.hf {
        RHF => println!("Number of Electrons:              {}", geometry.n_electrons),
        UHF => {
            println!(
                "Number of Electrons (Alpha):      {}",
                geometry.n_electrons_alpha
            );
            println!(
                "Number of Electrons (Beta):       {}",
                geometry.n_electrons_beta
            );
        }
        _ => panic!("Unsupported HF Type: {}", scf_input.hf),
    }
    println!("Basis Dimension:                  {}", basis.dim());

    println!("\nConverger\n---------");
    println!("Maximum Iterations:    {}", scf_input.max_iter);
    println!("Acceleration:          DIIS");

    println!("\nConvergence Thresholds\n----------------------");
    println!("Energy Change:         {:5.3e}", scf_input.e_threshold);
    println!("RMS:                   {:5.3e}", scf_input.rms_threshold);

    // read integrals from disk
    let h = FMatrix::retrieve(OneElectronKernel::HCore.to_filename());
    let s = FMatrix::retrieve(OneElectronKernel::Overlap.to_filename());

    // time
    let start_time = Instant::now();
    let eri = integrals.calc_two_electron_integral(TwoElectronKernel::ERI); // todo: verify
    println!("Time to calculate ERI: {:?}", start_time.elapsed());

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
