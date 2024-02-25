use serde::{Deserialize, Serialize};

use crate::geometry::atom::Atom;

use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
    str::FromStr,
};

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    atoms: Vec<Atom>,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>) -> Self {
        Self { atoms }
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn scale_coords(&mut self, factor: f64) {
        self.atoms
            .iter_mut()
            .for_each(|atom| atom.scale_coord(&factor));
    }

    /// Convert Molecule to xyz format
    fn to_xyz(&self) -> String {
        return format!("{}\n\n{}", self.atoms.len(), self);
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

impl FromStr for Molecule {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut lines = s.lines();

        // read number of atoms
        let n: usize = match lines.next() {
            Some(val) => val
                .trim()
                .parse()
                .map_err(|_| "Invalid number of coordinates")?,
            None => return Err("Missing xyz data"),
        };

        // skip comment line
        match lines.next() {
            Some(val) => (),
            None => return Err("This xyz string only contains one line"),
        };

        // read coordinates
        let mut atoms = Vec::<Atom>::with_capacity(n);
        for _ in 0..n {
            match lines.next() {
                Some(atom) => atoms.push(atom.trim().parse()?),
                None => {
                    return Err(
                        "The number of atoms does not agree with the number of coordinates found",
                    )
                }
            }
        }

        Ok(Self { atoms })
    }
}

impl Molecule {
    pub fn store(&self, name: &str) {
        let mut buffer =
            File::create(name.to_owned() + ".molecule").expect("Unable to create file");
        writeln!(buffer, "{}", self.to_xyz()).expect("Unable to write to file");
    }

    pub fn retrieve(name: &str) -> Self {
        let mut file =
            File::open(name.to_owned() + ".molecule").expect("Unable to open file for reading");
        let mut buffer = String::new();
        file.read_to_string(&mut buffer)
            .expect("Unable to read file");
        Self::from_str(&buffer).expect("Unable to deserialize a Molecule")
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
            Molecule::new(vec![
                Atom::new(Element::O, [0.1, 0.2, 3.0]),
                Atom::new(Element::H, [-3.3, -2.2, -1.1]),
            ],)
            .to_xyz(),
            r#"2

O     0.100000000    0.200000000    3.000000000
H    -3.300000000   -2.200000000   -1.100000000
"#
        );
    }

    #[test]
    fn deserialize_molecule() {
        use std::str::FromStr;

        // passing deserialisation
        assert_eq!(
            Molecule::from_str(
                r#"2

O     0.188972613    0.377945227    5.669178402
H    -6.236096242   -4.157397495   -2.078698747
"#
            ),
            Ok(Molecule {
                atoms: vec![
                    Atom::new(Element::O, [0.188972613, 0.377945227, 5.669178402]),
                    Atom::new(Element::H, [-6.236096242, -4.157397495, -2.078698747]),
                ],
            })
        );
        // failing deserialisations
        assert_eq!(
            Molecule::from_str("a"),
            Err("Invalid number of coordinates")
        );

        assert_eq!(Molecule::from_str(""), Err("Missing xyz data"));

        assert_eq!(
            Molecule::from_str(
                r#"2
"#
            ),
            Err("This xyz string only contains one line")
        );

        assert_eq!(
            Molecule::from_str(
                r#"2

O     0.188972613    0.377945227    5.669178402
"#
            ),
            Err("The number of atoms does not agree with the number of coordinates found")
        );
    }
}
