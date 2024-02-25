mod input;
use input::FerricInput;

use libferric::{
    gto_basis_sets::load_basis_set,
    gto_integrals::{
        integral_interface::IntegralInterface, one_electron::OneElectronKernel,
        two_electron::TwoElectronKernel,
    },
    misc::system::system_command,
};

use std::env::args;

fn print_banner() {
    println!(
        r#"
 _______  _______ .______      .______       __    ______
|   ____||   ____||   _  \     |   _  \     |  |  /      |
|  |__   |  |__   |  |_)  |    |  |_)  |    |  | |  ,----'
|   __|  |   __|  |      /     |      /     |  | |  |
|  |     |  |____ |  |\  \----.|  |\  \----.|  | |  `----.
|__|     |_______|| _| `._____|| _| `._____||__|  \______|"#
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

    // initialize basis set
    let basis = load_basis_set(&input.basis_set, input.geometry.molecule.atoms());
    basis.print_layout(input.geometry.molecule.atoms());
    basis.print_orca(input.geometry.molecule.atoms());

    // store stuff on disk
    input.geometry.molecule.store(&input.base_name);
    basis.store(&input.base_name);

    // todo: store "GBW" file on disk

    // calculate all necessary AO integrals for HF
    let integrals = IntegralInterface::new(&basis, input.geometry.molecule.atoms());
    integrals.calc_one_electron_integral(OneElectronKernel::Overlap);
    integrals.calc_one_electron_integral(OneElectronKernel::HCore);
    // integrals.calc_two_electron_integral(TwoElectronKernel::ERI); // todo: verify

    system_command("ferric_scf").expect("Something went wrong.");

    // remove all integrals
    integrals.remove();
}
