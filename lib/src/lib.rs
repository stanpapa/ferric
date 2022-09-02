pub mod basis;
pub mod constants;
mod elements;
pub mod geometry;
pub mod gto_basis_sets;
pub mod gto_integrals;
pub mod math;
pub mod misc;
pub mod orbitals;

extern crate blas;
extern crate blas_src;
extern crate lapack;
extern crate lapack_src;

use std::fmt::{Display, Formatter};
use std::str::FromStr;

// todo: move
pub enum HFType {
    RHF,
    UHF,
    ROHF,
    CASSCF,
}

impl Default for HFType {
    fn default() -> Self {
        HFType::RHF
    }
}

impl From<String> for HFType {
    fn from(s: String) -> Self {
        HFType::from_str(s.as_str()).unwrap()
    }
}

impl Display for HFType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            HFType::RHF => write!(f, "RHF"),
            HFType::UHF => write!(f, "UHF"),
            HFType::ROHF => write!(f, "ROHF"),
            HFType::CASSCF => write!(f, "CASSCF"),
        }
    }
}

impl FromStr for HFType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let hf = match s.to_uppercase().as_str() {
            "RHF" => HFType::RHF,
            "UHF" => HFType::UHF,
            "ROHF" => HFType::ROHF,
            "CASSCF" => HFType::CASSCF,
            _ => panic!("Unknown HF type {}", s),
        };

        Ok(hf)
    }
}
