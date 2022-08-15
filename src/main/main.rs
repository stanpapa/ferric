use ferric_lib::gto_bases;
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

    gto_bases::def2_tzvp::load_def2_tzvp();

    system_command("ferric-scf").expect("Something went wrong.");
}
