use std::str::FromStr;

pub mod driver;
mod hcore;

#[derive(Default)]
pub enum Guess {
    #[default]
    HCore,
}

impl FromStr for Guess {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "hcore" => Ok(Guess::HCore),
            _ => panic!("Unknown guess {}", s),
        }
    }
}
