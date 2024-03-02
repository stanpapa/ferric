use crate::{
    guess::{hcore, Guess},
    input::FerricInput,
};

pub fn driver(input: &FerricInput) {
    println!(
        r#"
▐▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▌
▐                                  ▌
▐      ____                        ▌
▐     / ___|_   _  ___  ___ ___    ▌
▐    | |  _| | | |/ _ \/ __/ __|   ▌
▐    | |_| | |_| |  __/\__ \__ \   ▌
▐     \____|\__,_|\___||___/___/   ▌
▐                                  ▌
▐                                  ▌
▐▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▌
"#
    );
    match input.guess {
        Guess::HCore => hcore::guess(&input.base_name, &input.scf.hf, &input.geometry),
    }
}
