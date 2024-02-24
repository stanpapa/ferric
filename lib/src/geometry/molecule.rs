use crate::{
    geometry::atom::Atom,
    linear_algebra::constants::{_ANG_AU, _AU_ANG, _BOHR_AU},
};

use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
};

use serde::{Deserialize, Serialize};
use toml::value::Table;

use super::Unit;

#[derive(Debug, Clone, Default, Deserialize, Serialize, PartialEq)]
pub struct Molecule {
    // basic information regarding a molecule (input)
    pub charge: i8,
    pub multiplicity: u8,

    #[serde(skip_serializing, skip_deserializing)]
    unit: Unit,

    // number of electrons (calculated)
    // #[serde(skip_serializing)]
    pub n_electrons: usize,
    // #[serde(skip_serializing)]
    pub n_electrons_alpha: usize,
    // #[serde(skip_serializing)]
    pub n_electrons_beta: usize,

    atoms: Vec<Atom>,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>, charge: i8, multiplicity: u8) -> Molecule {
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

        let mut mol = Molecule {
            charge,
            multiplicity,
            unit: Unit::Ångström,
            atoms,
            n_electrons,
            n_electrons_alpha,
            n_electrons_beta,
        };

        // convert coordinates to atomic units
        mol.scale_coords(Unit::AtomicUnits);

        mol
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn scale_coords(&mut self, unit: Unit) {
        for atom in &mut self.atoms {
            match unit {
                Unit::Ångström => match self.unit {
                    Unit::Ångström => (),
                    Unit::AtomicUnits => atom.scale_coords(&_AU_ANG),
                    Unit::Bohr => todo!(),
                },
                Unit::AtomicUnits => match self.unit {
                    Unit::Ångström => atom.scale_coords(&_ANG_AU),
                    Unit::AtomicUnits => (),
                    Unit::Bohr => todo!(),
                },
                Unit::Bohr => match self.unit {
                    Unit::Ångström => todo!(),
                    Unit::AtomicUnits => todo!(),
                    Unit::Bohr => (),
                },
            }
        }
        self.unit = unit;
    }

    pub fn print(&self, unit: Unit) {
        // ---------------------------------
        // CARTESIAN COORDINATES (ANGSTROEM)
        // ---------------------------------
        //   O      0.000000    0.000000   -0.119015
        //   H      0.768550    0.000000    0.476060
        //   H     -0.768550    0.000000    0.476060

        // ----------------------------
        // CARTESIAN COORDINATES (A.U.)
        // ----------------------------
        //   NO LB      ZA    FRAG     MASS         X           Y           Z
        //    0 O     8.0000    0    15.999    0.000000    0.000000   -0.224906
        //    1 H     1.0000    0     1.008    1.452350    0.000000    0.899624
        //    2 H     1.0000    0     1.008   -1.452350    0.000000    0.899624
        let mut mol = self.clone();
        match unit {
            Unit::Ångström => {
                println!(
                    r#"
--------------------------------
CARTESIAN COORDINATES (ÅNGSTRÖM)
--------------------------------
"#
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
        mol.scale_coords(unit);
        println!("{}", mol);
    }
}

impl Display for Molecule {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for atom in &self.atoms {
            write!(f, "{}", atom)?;
        }

        Ok(())
    }
}

impl Molecule {
    pub fn from_table(table: Table) -> Self {
        let mut charge = 0;
        let mut mult = 1;
        let mut coordinates: Vec<Atom> = vec![];

        for (key, val) in table.iter() {
            match key.to_lowercase().as_str() {
                "charge" => charge = val.as_integer().unwrap().try_into().unwrap(),
                "multiplicity" => mult = val.as_integer().unwrap().try_into().unwrap(),
                "coordinates" => match val.as_array() {
                    Some(array) => {
                        coordinates = array.iter().map(|i| Atom::from_array(i)).collect()
                    }
                    None => panic!("Coordinates input is not valid"),
                },
                _ => panic!("Unknown option: {}", key),
            }
        }

        Self::new(coordinates, charge, mult)
    }

    pub fn store(&self, name: &str) {
        let mut buffer =
            File::create(name.to_owned() + ".molecule").expect("Unable to create file");
        write!(
            buffer,
            "{}",
            toml::to_string(self).expect("Unable to serialize Molecule")
        )
        .expect("Unable to write to file");
    }

    pub fn retrieve(name: &str) -> Self {
        let mut file =
            File::open(name.to_owned() + ".molecule").expect("Unable to open file for reading");
        let mut buffer = String::new();
        file.read_to_string(&mut buffer)
            .expect("Unable to read file");
        toml::from_str(&buffer).expect("Unable to deserialize a Molecule")
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::molecule::Atom;
    use crate::geometry::molecule::Molecule;
    use crate::misc::elements::Element;

    #[test]
    fn serialize_molecule() {
        assert_eq!(
            toml::to_string(&Molecule::new(
                vec![
                    Atom::new(Element::O, [0.1, 0.2, 3.0]),
                    Atom::new(Element::H, [-3.3, -2.2, -1.1]),
                ],
                -1,
                1
            ))
            .unwrap(),
            r#"charge = -1
multiplicity = 1

[[atoms]]
el = "O"
origin = [0.18897261339212518, 0.37794522678425035, 5.669178401763755]

[[atoms]]
el = "H"
origin = [-6.23609624194013, -4.157397494626754, -2.078698747313377]
"#
        );
    }

    #[test]
    fn deserialize_molecule() {
        assert_eq!(
            toml::from_str::<Molecule>(
                r#"charge = -1
multiplicity = 1

[[atoms]]
el = "O"
origin = [0.18897261339212518, 0.37794522678425035, 5.669178401763755]

[[atoms]]
el = "H"
origin = [-6.23609624194013, -4.157397494626754, -2.078698747313377]
"#
            )
            .unwrap(),
            Molecule::new(
                vec![
                    Atom::new(Element::O, [0.1, 0.2, 3.0]),
                    Atom::new(Element::H, [-3.3, -2.2, -1.1]),
                ],
                -1,
                1
            )
        );
    }
}
