pub mod data;
pub mod geometry;
pub mod gto_basis_sets;
pub mod gto_integrals;
pub mod linear_algebra;
pub mod misc;

use serde::{Deserialize, Serialize};

extern crate blas_src;
extern crate cblas;
extern crate lapack_src;
extern crate lapacke;

use std::{
    fmt::{Display, Formatter},
    str::FromStr,
};

// todo: move
#[derive(PartialEq, Serialize, Deserialize)]
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

#[cfg(test)]
mod tests {
    use crate::HFType;

    #[test]
    fn serialize() {
        assert_eq!(toml::to_string(&HFType::RHF).unwrap(), r#""RHF""#);
        assert_eq!(toml::to_string(&HFType::UHF).unwrap(), r#""UHF""#);
        assert_eq!(toml::to_string(&HFType::ROHF).unwrap(), r#""ROHF""#);
        assert_eq!(toml::to_string(&HFType::CASSCF).unwrap(), r#""CASSCF""#);
    }

    #[test]
    fn deserialize() {
        assert_eq!(r#""RHF""#, toml::to_string(&HFType::RHF).unwrap());
        assert_eq!(r#""UHF""#, toml::to_string(&HFType::UHF).unwrap());
        assert_eq!(r#""ROHF""#, toml::to_string(&HFType::ROHF).unwrap());
        assert_eq!(r#""CASSCF""#, toml::to_string(&HFType::CASSCF).unwrap(),);
    }
}
