mod input;

use crate::input::FerricInput;
use ferric_lib::gto_basis_sets::load_basis_set;
use ferric_lib::gto_integrals::integral_interface::IntegralInterface;
use ferric_lib::gto_integrals::nuclear_repulsion::nuclear_repulsion;
use ferric_lib::gto_integrals::one_electron::OneElectronKernel;
use ferric_lib::math::matrix::FMatrix;
use ferric_lib::math::traits::Diagonalize;
use ferric_lib::misc::system::system_command;
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

    let input = FerricInput::new(&mut args());

    let basis = load_basis_set(input.basis_set(), input.molecule().atoms());

    let integrals = IntegralInterface::new(&basis, input.molecule().atoms());

    let nuclear_repulsion = nuclear_repulsion(input.molecule().atoms());
    println!("Vnn: {}", nuclear_repulsion);

    let overlap = integrals.calc_one_electron_integral(OneElectronKernel::Overlap);
    println!("Overlap:\n{}", overlap); // verified aginst EasyIntegrals

    let kinetic = integrals.calc_one_electron_integral(OneElectronKernel::Kinetic);
    println!("Kinetic\n{}", kinetic); // verified aginst EasyIntegrals

    let v = integrals.calc_one_electron_integral(OneElectronKernel::NuclearAttraction);
    println!("Nuclear attraction\n{}", v); // verified aginst EasyIntegrals

    let h = kinetic + v;
    println!("One-electron matrix:\n{}", h);

    // // diagonalize overlap
    // let (eigenvalues, eigenvectors) = overlap.diagonalize();
    // let mut lambda_12 = FMatrix::zero(eigenvalues.size(), eigenvalues.size());
    // for i in 0..lambda_12.cols() {
    //     lambda_12[(i, i)] = 1.0 / eigenvalues[i].sqrt();
    // }
    //
    // let s12 = eigenvectors.clone() * (lambda_12 * eigenvectors.transposed());
    // println!("S^(-1/2):\n{}", s12); // verified aginst EasyIntegrals
    //
    // // core fock matrix: F0' = S^{-1/2}^T H S^{-1/2}
    // // s-12 = symmetric
    // let h_mat = FMatrix::from(h);
    // let f0 = s12.clone() * (h_mat * s12.clone());
    // println!("F0:\n{}", f0); // verified aginst EasyIntegrals
    //
    // // diagonlize F0': C0'^T F0' C0' = eps
    // let (_eps, c0prime) = f0.diagonalize();
    //
    // let c0 = s12 * c0prime;
    // println!("C0:\n{}", c0); // verified aginst EasyIntegrals

    system_command("ferric_scf").expect("Something went wrong.");
}
