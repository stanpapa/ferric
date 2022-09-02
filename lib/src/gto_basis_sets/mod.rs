use crate::basis::Basis;
use crate::geometry::Atom;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

mod def2_tzvp;
mod sto_3g;

use def2_tzvp::load_def2_tzvp;
use sto_3g::load_sto_3g;

#[allow(non_camel_case_types)]
#[derive(Debug, PartialEq)]
pub enum BasisSet {
    sto_3g,
    def2_tzvp,
}

impl Default for BasisSet {
    fn default() -> Self {
        BasisSet::sto_3g
    }
}

impl Display for BasisSet {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s;
        match self {
            BasisSet::sto_3g => s = "sto-3g",
            BasisSet::def2_tzvp => s = "def2-tzvp",
        }
        write!(f, "{}", s)
    }
}

pub fn load_basis_set(basis_set: &BasisSet, atoms: &[Atom]) -> Basis {
    match basis_set {
        BasisSet::sto_3g => load_sto_3g(atoms),
        BasisSet::def2_tzvp => load_def2_tzvp(atoms),
    }
}

impl From<String> for BasisSet {
    fn from(s: String) -> Self {
        BasisSet::from_str(s.as_str()).unwrap()
    }
}

impl FromStr for BasisSet {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let basis = match s.to_lowercase().as_str() {
            "sto-3g" => BasisSet::sto_3g,
            "def2-tzvp" => BasisSet::def2_tzvp,
            _ => panic!("Unknown basis set {}", s),
        };

        Ok(basis)
    }
}

#[cfg(test)]
mod tests {
    use crate::gto_basis_sets::BasisSet;
    use std::str::FromStr;

    #[test]
    fn from() {
        assert_eq!(BasisSet::from("sto-3g".to_string()), BasisSet::sto_3g);
        assert_eq!(BasisSet::from("sTO-3g".to_string()), BasisSet::sto_3g);
        assert_eq!(BasisSet::from("STO-3G".to_string()), BasisSet::sto_3g);
    }

    #[test]
    fn from_str() {
        assert_eq!(BasisSet::from_str("sto-3g"), Ok(BasisSet::sto_3g));
        assert_eq!(BasisSet::from_str("sTO-3g"), Ok(BasisSet::sto_3g));
        assert_eq!(BasisSet::from_str("STO-3G"), Ok(BasisSet::sto_3g));
    }

    #[test]
    #[should_panic]
    #[allow(unused_must_use)]
    fn invalid_from_str() {
        BasisSet::from_str("sto_3g");
    }

    #[test]
    #[should_panic]
    #[allow(unused_must_use)]
    fn invalid_from() {
        BasisSet::from("sto_3g".to_string());
    }
}
