use libferric::{
    geometry::{Geometry, Unit},
    gto_basis_sets::BasisSet,
};

use crate::{guess::Guess, scf::input::SCFInput};

use std::{env::Args, fs, str::FromStr};

use serde_yaml::Value;

#[derive(Default)]
pub struct FerricInput {
    pub base_name: String,

    // Basis
    pub basis_set: BasisSet,

    // Geometry
    pub geometry: Geometry,

    // Modules
    pub guess: Guess,
    pub scf: SCFInput,
}

// todo: use Serialize/Deserialize?
impl FerricInput {
    pub fn new(input: &mut Args) -> Self {
        let mut ferric_input = Self::default();

        // determine base name of file
        let input_file = input.nth(1).expect("No input file given");
        ferric_input.base_name = input_file
            .strip_suffix(".yaml")
            .expect("Input file does not end with .yaml")
            .to_string();

        // read input file as a toml::Table
        let contents = fs::read_to_string(input_file).expect("Not able to read input file");
        let inp: Value = serde_yaml::from_str(&contents).unwrap();
        for (key, value) in inp.as_mapping().unwrap() {
            match key.as_str().unwrap().to_lowercase().as_str() {
                "guess" => ferric_input.parse_guess(value),
                "scf" => ferric_input.scf = SCFInput::parse(value),
                "basis" => ferric_input.parse_basis(value),
                "geometry" => ferric_input.parse_geometry(value),
                _ => panic!("Invalid block {:?}", key),
            }
        }

        ferric_input
    }
}

impl FerricInput {
    fn parse_guess(&mut self, input: &Value) {
        match input {
            Value::String(s) => self.guess = Guess::from_str(s).unwrap(),
            _ => panic!("Invalid guess {:?}", input),
        }
    }

    fn parse_basis(&mut self, input: &Value) {
        match input {
            Value::String(s) => self.basis_set = BasisSet::from_str(s).unwrap(),
            _ => panic!("Invalid basis {:?}", input),
        }
    }

    pub fn parse_geometry(&mut self, value: &Value) {
        // todo!("Error management")
        match value.as_mapping() {
            Some(mapping) => {
                let charge = match mapping.get("charge") {
                    Some(charge) => charge.as_i64().unwrap(),
                    None => panic!("No charge found"),
                };
                let multiplicity = match mapping.get("multiplicity") {
                    Some(multiplicity) => multiplicity.as_u64().unwrap(),
                    None => panic!("No multiplicity found"),
                };
                let unit = match mapping.get("unit") {
                    Some(unit) => Unit::from_str(unit.as_str().unwrap()).unwrap(),
                    None => Unit::Ångström,
                };
                let atoms = match mapping.get("xyz") {
                    Some(atoms) => atoms.as_str().unwrap(),
                    None => panic!("No coordinates found"),
                };
                self.geometry = Geometry::new(
                    atoms.lines().map(|atom| atom.parse().unwrap()).collect(),
                    charge as i8,
                    multiplicity as u8,
                    unit,
                );
            }
            None => panic!("Invalid geometry input {:?}", value),
        }
    }
}
