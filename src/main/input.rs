use ferric_lib::geometry::Molecule;
use ferric_lib::gto_basis_sets::BasisSet;
use ferric_lib::HFType;
use std::env::Args;
use std::fs;
use std::str::FromStr;

use toml::value::Table;
use toml::Value;

#[derive(Default)]
pub struct FerricInput {
    base_name: String,

    hf: HFType,
    basis_set: BasisSet,

    // Geometry
    mol: Molecule,
}

// todo: use Serialize/Deserialize?
impl FerricInput {
    pub fn new(input: &mut Args) -> Self {
        let mut ferric_input = Self::default();

        // determine base name of file
        let input_file = input.nth(1).expect("No input file given");
        ferric_input.base_name = input_file
            .strip_suffix(".toml")
            .expect("Input file does not end with .toml")
            .to_string();

        // read input file as a toml::Table
        let contents = fs::read_to_string(input_file).expect("Not able to read input file");
        let inp = contents.parse::<Value>().unwrap();
        let table = match inp.as_table() {
            Some(t) => t,
            None => panic!("Not able to parse input"),
        };

        // parse subsections of input
        for (key, value) in table.iter() {
            match key.to_lowercase().as_str() {
                "main" => ferric_input.parse_main(value),
                "geometry" => ferric_input.parse_geometry(value),
                _ => panic!("Unknown block {}", key),
            }
        }

        ferric_input
    }

    // pub fn base_name(&self) -> &str {
    //     &self.base_name
    // }

    pub fn basis_set(&self) -> &BasisSet {
        &self.basis_set
    }

    pub fn molecule(&self) -> &Molecule {
        &self.mol
    }
}

impl FerricInput {
    fn check_if_table(value: &Value) -> Table {
        value.as_table().expect("Not able to parse input").clone()
    }

    pub fn parse_main(&mut self, value: &Value) {
        let table = Self::check_if_table(value);

        for (key, val) in table.iter() {
            match key.to_lowercase().as_str() {
                "hf" => self.hf = HFType::from_str(val.as_str().unwrap()).unwrap(),
                "basis" => self.basis_set = BasisSet::from_str(val.as_str().unwrap()).unwrap(),
                _ => panic!("Unknown option: {}", key),
            }
        }
    }

    pub fn parse_geometry(&mut self, value: &Value) {
        let table = Self::check_if_table(value);
        self.mol = Molecule::from_table(table);
    }
}
