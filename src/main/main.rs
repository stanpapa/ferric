use ferric_lib::geometry::Molecule;
use ferric_lib::gto_bases::sto_3g::load_sto_3g;
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

    load_sto_3g(mol.atoms());

    system_command("ferric_scf").expect("Something went wrong.");
}
