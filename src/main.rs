// main directory
mod input;

// submodules
mod guess;
mod scf;

use input::FerricInput;

use libferric::{
    gto_basis_sets::load_basis_set,
    gto_integrals::{
        integral_interface::IntegralInterface, one_electron::OneElectronKernel,
        two_electron::TwoElectronKernel,
    },
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

    // --------------------------------------------------
    // read input file
    // --------------------------------------------------
    let input = FerricInput::new(&mut args());

    // --------------------------------------------------
    // print geometry
    // --------------------------------------------------
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
    println!(
        r#"
▐▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▌
▐                                               ▌
▐    ___       _                       _        ▌
▐   |_ _|_ __ | |_ ___  __ _ _ __ __ _| |___    ▌
▐    | || '_ \| __/ _ \/ _` | '__/ _` | / __|   ▌
▐    | || | | | ||  __/ (_| | | | (_| | \__ \   ▌
▐   |___|_| |_|\__\___|\__, |_|  \__,_|_|___/   ▌
▐                      |___/                    ▌
▐                                               ▌
▐▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▌
"#
    );
    let integrals = IntegralInterface::new(&basis, input.geometry.molecule.atoms());
    let _h = integrals.calc_one_electron_integral(OneElectronKernel::Overlap);
    let _s = integrals.calc_one_electron_integral(OneElectronKernel::HCore);
    let _eri = integrals.calc_two_electron_integral(TwoElectronKernel::ERI);

    // --------------------------------------------------
    // Guess
    // --------------------------------------------------
    guess::driver::driver(&input).expect("Guess could be constructed succesfully");

    // --------------------------------------------------
    // SCF Calculation
    // --------------------------------------------------
    scf::driver::driver(&input.base_name, input.scf)
        .expect("SCF calculation did not finish succesfully");

    // --------------------------------------------------
    // clean-up
    // --------------------------------------------------
    integrals.remove();
    // todo!() remove basis/geometry
}
