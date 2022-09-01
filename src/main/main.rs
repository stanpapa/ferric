use ferric_lib::geometry::Molecule;
use ferric_lib::gto_basis_sets::sto_3g::load_sto_3g;
use ferric_lib::gto_integrals::integral_interface::IntegralInterface;
use ferric_lib::gto_integrals::nuclear_repulsion::nuclear_repulsion;
use ferric_lib::gto_integrals::one_electron::OneElectronKernel;
use ferric_lib::misc::system::system_command;

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

    let mol = Molecule::default();

    let basis = load_sto_3g(mol.atoms());

    let integrals = IntegralInterface::new(&basis, mol.atoms());

    let nuclear_repulsion = nuclear_repulsion(&mol);
    println!("Vnn: {}", nuclear_repulsion);

    let overlap = integrals.calc_one_electron_integral(OneElectronKernel::Overlap);
    println!("Overlap:\n{}", overlap);

    let kinetic = integrals.calc_one_electron_integral(OneElectronKernel::Kinetic);
    println!("Kinetic\n{}", kinetic);

    let v = integrals.calc_one_electron_integral(OneElectronKernel::NuclearAttracton);
    println!("Nuclear attraction\n{}", v);

    system_command("ferric_scf").expect("Something went wrong.");
}
