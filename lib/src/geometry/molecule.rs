use crate::{geometry::atom::Atom, linear_algebra::constants::_ANG_AU};

use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
};

use serde::{Deserialize, Serialize};
use toml::value::Table;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Molecule {
    pub charge: i8,
    pub multiplicity: u8,
    atoms: Vec<Atom>,
    // num_electrons_alpha: usize,
    // num_electrons_beta: usize,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>, charge: i8, multiplicity: u8) -> Molecule {
        // let mut num_electrons: i8 = -charge;
        // num_electrons += atoms.iter().map(|a| a.z()).sum::<u8>() as i8;
        //
        // if num_electrons <= 0 {
        //     panic!("Invalid charge of {}.", charge);
        // }
        //
        // // Mult = 2Ms + 1 thus the number of unpaired electrons (taken as α) is Mult-1 = 2Ms
        // let alpha_excess = multiplicity as usize - 1;
        //
        // // the number of β electrons must be equal the number of doubly occupied orbitals
        // // Ndo = (nelec - n_of_unpaired_elec)/2 this must be an integer
        // // if isodd(num_electrons) != isodd(alpha_excess) {
        // if (num_electrons % 2 == 0) && (alpha_excess % 2 == 1) {
        //     panic!(
        //         "Incompatible charge {} and multiplicity {}",
        //         charge, multiplicity
        //     );
        // }
        //
        // let num_electrons_beta = (num_electrons as usize - alpha_excess) / 2;
        // let num_electrons_alpha = num_electrons as usize - num_electrons_beta;

        let mut mol = Molecule {
            charge,
            multiplicity,
            atoms,
            // num_electrons_alpha,
            // num_electrons_beta,
        };

        mol.scale_coords(&_ANG_AU);

        mol
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn scale_coords(&mut self, factor: &f64) {
        for atom in &mut self.atoms {
            atom.scale_coords(factor);
        }
    }
}

// of course, the default molecule has to be H2O
impl Default for Molecule {
    fn default() -> Self {
        Self::new(vec![], 0, 0)
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
