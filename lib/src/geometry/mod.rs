use std::{collections::HashMap, fmt::Display, str::FromStr};

use crate::linear_algebra::constants::{ANG_AU, ANG_BOHR, AU_ANG, AU_BOHR, BOHR_ANG, BOHR_AU};

use self::{atom::Atom, molecule::Molecule};

pub mod atom;
pub mod molecule;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Hash)]
pub enum Unit {
    #[default]
    Ångström,
    AtomicUnits,
    Bohr,
}

impl Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unit::Ångström => write!(f, "ångström")?,
            Unit::AtomicUnits => write!(f, "au")?,
            Unit::Bohr => write!(f, "bohr")?,
        }

        Ok(())
    }
}

impl FromStr for Unit {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str().trim() {
            "ångström" | "angstroem" | "" => Ok(Unit::Ångström),
            "au" => Ok(Unit::AtomicUnits),
            "bohr" => Ok(Unit::Bohr),
            _ => Err("Unit::from_str: Invalid unit"),
        }
    }
}

#[derive(Default)]
pub struct Geometry {
    pub molecule: Molecule,
    unit: Unit,

    pub charge: i8,
    pub multiplicity: u8,

    // number of electrons (calculated)
    pub n_electrons: usize,
    pub n_electrons_alpha: usize,
    pub n_electrons_beta: usize,
}

impl Geometry {
    pub fn new(atoms: Vec<Atom>, charge: i8, multiplicity: u8, unit: Unit) -> Self {
        // calculate number of electrons
        let n_electrons: isize =
            atoms.iter().map(|a| isize::from(a.z())).sum::<isize>() - isize::from(charge);

        if n_electrons <= 0 {
            panic!("Invalid charge of {}. No electrons present", charge);
        }

        // Mult = 2Ms + 1 thus the number of unpaired electrons (taken as α) is Mult-1 = 2Ms
        let alpha_excess = usize::from(multiplicity) - 1;

        // the number of β electrons must be equal the number of doubly occupied orbitals
        // Ndo = (nelec - n_of_unpaired_elec)/2 this must be an integer
        if (n_electrons % 2 == 1) && (multiplicity % 2 == 1) {
            panic!(
                "Both multiplicity ({}) and number of electrons ({}) are odd, which is impossible",
                multiplicity, n_electrons
            );
        }
        let n_electrons = n_electrons as usize;

        let n_electrons_beta = (n_electrons - alpha_excess) / 2;
        let n_electrons_alpha = n_electrons - n_electrons_beta;

        let mut geometry = Self {
            molecule: Molecule::new(atoms),
            unit,
            charge,
            multiplicity,
            n_electrons,
            n_electrons_alpha,
            n_electrons_beta,
        };

        // convert coordinates to atomic units
        geometry.convert_unit(Unit::AtomicUnits);

        geometry
    }

    pub fn convert_unit(&mut self, unit: Unit) {
        if self.unit == unit {
            return;
        }

        let convert: HashMap<(Unit, Unit), f64> = [
            ((Unit::Ångström, Unit::AtomicUnits), ANG_AU),
            ((Unit::Ångström, Unit::Bohr), ANG_BOHR),
            ((Unit::AtomicUnits, Unit::Ångström), AU_ANG),
            ((Unit::AtomicUnits, Unit::Bohr), AU_BOHR),
            ((Unit::Bohr, Unit::Ångström), BOHR_ANG),
            ((Unit::Bohr, Unit::AtomicUnits), BOHR_AU),
        ]
        .into_iter()
        .collect();

        self.molecule
            .scale_coords(*convert.get(&(self.unit, unit)).unwrap());

        self.unit = unit;
    }

    pub fn print_coords(&self, unit: Unit) {
        let mut tmp = self.molecule.clone();
        match unit {
            Unit::Ångström => {
                println!(
                    r#"
--------------------------------
CARTESIAN COORDINATES (ANGSTROM)
--------------------------------
"#,
                );
            }
            Unit::Bohr => {
                println!(
                    r#"
----------------------------
CARTESIAN COORDINATES (BOHR)
----------------------------
"#
                );
            }
            Unit::AtomicUnits => {
                println!(
                    r#"
----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
"#
                );
            }
        }

        let convert: HashMap<(Unit, Unit), f64> = [
            ((Unit::Ångström, Unit::Ångström), 1.0),
            ((Unit::Ångström, Unit::AtomicUnits), ANG_AU),
            ((Unit::Ångström, Unit::Bohr), ANG_BOHR),
            ((Unit::AtomicUnits, Unit::Ångström), AU_ANG),
            ((Unit::AtomicUnits, Unit::AtomicUnits), 1.0),
            ((Unit::AtomicUnits, Unit::Bohr), AU_BOHR),
            ((Unit::Bohr, Unit::Ångström), BOHR_ANG),
            ((Unit::Bohr, Unit::AtomicUnits), BOHR_AU),
            ((Unit::Bohr, Unit::Bohr), 1.0),
        ]
        .into_iter()
        .collect();

        tmp.scale_coords(*convert.get(&(self.unit, unit)).unwrap());

        println!("{}", tmp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_to_string() {
        assert_eq!(Unit::Ångström.to_string(), String::from("ångström"));
        assert_eq!(Unit::AtomicUnits.to_string(), String::from("au"));
        assert_eq!(Unit::Bohr.to_string(), String::from("bohr"));
    }

    #[test]
    fn unit_from_str() {
        assert_eq!(Unit::from_str("ångström "), Ok(Unit::Ångström));
        assert_eq!(Unit::from_str("aNgstroem"), Ok(Unit::Ångström));
        assert_eq!(Unit::from_str(""), Ok(Unit::Ångström));
        assert_eq!(Unit::from_str("  "), Ok(Unit::Ångström));
        assert_eq!(Unit::from_str(" au"), Ok(Unit::AtomicUnits));
        assert_eq!(Unit::from_str(" bohr "), Ok(Unit::Bohr));
    }
}
