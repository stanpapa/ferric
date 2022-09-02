use crate::constants::_ANG_AU;
use crate::elements::Element;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

use toml::value::{Table, Value};

pub struct Molecule {
    atoms: Vec<Atom>,
    charge: i8,
    multiplicity: u8,
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
            atoms,
            charge,
            multiplicity,
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

#[derive(Clone)]
pub struct Atom {
    el: Element,
    origin: [f64; 3],
}

impl Atom {
    pub fn new(el: Element, origin: [f64; 3]) -> Self {
        Self { el, origin }
    }

    pub fn z(&self) -> u8 {
        self.el.atomic_number()
    }

    pub fn mass(&self) -> f32 {
        self.el.mass()
    }

    pub fn origin(&self) -> &[f64; 3] {
        &self.origin
    }
    //
    // pub fn origin_mut(&mut self) -> &mut [f64; 3] {
    //     &mut self.origin
    // }

    pub fn scale_coords(&mut self, factor: &f64) {
        for c in &mut self.origin {
            *c *= factor;
        }
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}: {:14.9}  {:14.9} {:14.9}",
            self.el, self.origin[0], self.origin[1], self.origin[2]
        )
    }
}

impl Atom {
    fn from_array(value: &Value) -> Self {
        match value.as_array() {
            Some(array) => Atom {
                el: Element::from_str(array[0].as_str().unwrap()).unwrap(),
                origin: array[1..=3]
                    .iter()
                    .map(|v| v.as_float().expect("coordinate is not a float"))
                    .collect::<Vec<f64>>()
                    .try_into()
                    .unwrap(),
            },
            None => panic!("Invalid atom array"),
        }
    }
}

pub fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];

    (dx * dx + dy * dy + dz * dz).sqrt()
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
}
