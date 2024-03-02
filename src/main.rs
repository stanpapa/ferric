// main directory
mod input;

// submodules
mod guess;
mod scf;

use input::FerricInput;

use libferric::{
    gto_basis_sets::load_basis_set,
    gto_integrals::{integral_interface::IntegralInterface, one_electron::OneElectronKernel},
};

use std::env::args;

fn print_banner() {
    println!(
        r#"
 _____                                    _____ 
( ___ )----------------------------------( ___ )
 |   |                                    |   | 
 |   |                                    |   | 
 |   |      _____              _          |   | 
 |   |     |  ___|__ _ __ _ __(_) ___     |   | 
 |   |     | |_ / _ \ '__| '__| |/ __|    |   | 
 |   |     |  _|  __/ |  | |  | | (__     |   | 
 |   |     |_|  \___|_|  |_|  |_|\___|    |   | 
 |   |                                    |   | 
 |   |                                    |   | 
 |___|                                    |___| 
(_____)----------------------------------(_____)
"#
    );
}

fn main() {
    print_banner();

    // read input file
    let input = FerricInput::new(&mut args());

    input
        .geometry
        .print_coords(libferric::geometry::Unit::Ångström);
    input
        .geometry
        .print_coords(libferric::geometry::Unit::AtomicUnits);
    input.geometry.store(&input.base_name);

    // --------------------------------------------------
    // initialize basis set
    // --------------------------------------------------
    let basis = load_basis_set(&input.basis_set, input.geometry.molecule.atoms());
    basis.print_layout(input.geometry.molecule.atoms());
    basis.print_orca(input.geometry.molecule.atoms());
    basis.store(&input.base_name);

    // todo: store "GBW" file on disk

    // --------------------------------------------------
    // calculate all necessary AO integrals for HF
    // --------------------------------------------------
    let integrals = IntegralInterface::new(&basis, input.geometry.molecule.atoms());
    integrals.calc_one_electron_integral(OneElectronKernel::Overlap);
    integrals.calc_one_electron_integral(OneElectronKernel::HCore);

    // --------------------------------------------------
    // Guess
    // --------------------------------------------------
    guess::driver::driver(&input);

    // --------------------------------------------------
    // SCF Calculation
    // --------------------------------------------------
    scf::driver::driver(&input.base_name, input.scf)
        .expect("SCF calculation did not finish succesfully");

    // remove all integrals
    integrals.remove();
    // todo!() remove basis/geometry
}
