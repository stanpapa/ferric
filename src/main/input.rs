use libferric::geometry::molecule::Molecule;
use libferric::geometry::{Geometry, Unit};
use libferric::gto_basis_sets::BasisSet;
use libferric::HFType;
use std::env::Args;
use std::fs;
use std::str::FromStr;

use serde_yaml::Value;

#[derive(Default)]
pub struct FerricInput {
    pub base_name: String,

    pub hf: HFType,
    pub basis_set: BasisSet,

    // Geometry
    pub geometry: Geometry,
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
                "scf" => ferric_input.parse_scf(value),
                "basis" => ferric_input.parse_basis(value),
                "geometry" => ferric_input.parse_geometry(value),
                _ => panic!("Invalid block {:?}", key),
            }
        }

        ferric_input
    }
}

impl FerricInput {
    fn parse_scf(&mut self, input: &Value) {
        for (key, value) in input.as_mapping().unwrap() {
            match key.as_str().unwrap().to_lowercase().as_str() {
                "hf" => self.hf = HFType::from_str(value.as_str().unwrap()).unwrap(),
                _ => panic!("Unknown option: {:?}", key),
            }
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
