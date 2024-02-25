pub mod data;
pub mod geometry;
pub mod gto_basis_sets;
pub mod gto_integrals;
pub mod linear_algebra;
pub mod misc;

extern crate blas_src;
extern crate cblas;
extern crate lapack_src;
extern crate lapacke;

use std::{
    fmt::{Display, Formatter},
    str::FromStr,
};

// todo: move
#[derive(Debug, Default, PartialEq)]
pub enum HFType {
    #[default]
    RHF,
    UHF,
    ROHF,
    CASSCF,
}

impl From<String> for HFType {
    fn from(s: String) -> Self {
        HFType::from_str(s.as_str()).unwrap()
    }
}

impl FromStr for HFType {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let hf = match s.to_uppercase().as_str() {
            "RHF" => HFType::RHF,
            "UHF" => HFType::UHF,
            "ROHF" => HFType::ROHF,
            "CASSCF" => HFType::CASSCF,
            _ => return Err("Unknown HF type"),
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
    use super::HFType;

    #[test]
    fn serialize() {
        assert_eq!(HFType::RHF.to_string(), "RHF");
        assert_eq!(HFType::UHF.to_string(), "UHF");
        assert_eq!(HFType::ROHF.to_string(), "ROHF");
        assert_eq!(HFType::CASSCF.to_string(), "CASSCF");
    }

    #[test]
    fn deserialize() {
        use std::str::FromStr;

        assert_eq!(HFType::from_str("RHF"), Ok(HFType::RHF));
        assert_eq!(HFType::from_str("UHF"), Ok(HFType::UHF));
        assert_eq!(HFType::from_str("ROHF"), Ok(HFType::ROHF));
        assert_eq!(HFType::from_str("CASSCF"), Ok(HFType::CASSCF));
    }
}
