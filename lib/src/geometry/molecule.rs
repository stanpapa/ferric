use crate::{
    geometry::atom::Atom,
    linear_algebra::constants::{_ANG_AU, _AU_ANG, _BOHR_AU},
};

use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
    str::FromStr,
};

use super::Unit;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    // basic information regarding a molecule (input)
    pub charge: i8,
    pub multiplicity: u8,

    unit: Unit,

    // number of electrons (calculated)
    pub n_electrons: usize,
    pub n_electrons_alpha: usize,
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

    /// Convert Molecule to xyz format
    fn to_xyz(&self) -> String {
        return format!("{}\n{}\n{}", self.atoms.len(), self.unit, self);
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

        let unit: Unit = match lines.next() {
            Some(val) => val.trim().parse().map_err(|_| "Invalid unit")?,
            None => return Err("This xyz string only contains one line"),
        };

        let mut atoms = Vec::<Atom>::with_capacity(n);
        for _ in 0..n {
            // atoms.push(lines.next().unwrap().trim().parse()?);
            match lines.next() {
                Some(atom) => atoms.push(atom.trim().parse()?),
                None => {
                    return Err(
                        "The number of atoms does not agree with the number of coordinates found",
                    )
                }
            }
        }

        Ok(Self {
            charge: i8::MAX,
            multiplicity: u8::MAX,
            unit,
            n_electrons: 0,
            n_electrons_alpha: 0,
            n_electrons_beta: 0,
            atoms,
        })
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
            Molecule::new(
                vec![
                    Atom::new(Element::O, [0.1, 0.2, 3.0]),
                    Atom::new(Element::H, [-3.3, -2.2, -1.1]),
                ],
                -1,
                1
            )
            .to_xyz(),
            r#"2
au
O     0.188972613    0.377945227    5.669178402
H    -6.236096242   -4.157397495   -2.078698747
"#
        );
    }

    #[test]
    fn deserialize_molecule() {
        use super::Unit;
        use std::str::FromStr;

        // passing deserialisation
        assert_eq!(
            Molecule::from_str(
                r#"2
au
O     0.188972613    0.377945227    5.669178402
H    -6.236096242   -4.157397495   -2.078698747
"#
            ),
            Ok(Molecule {
                charge: i8::MAX,
                multiplicity: u8::MAX,
                unit: Unit::AtomicUnits,
                n_electrons: 0,
                n_electrons_alpha: 0,
                n_electrons_beta: 0,
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
foo
O     0.188972613    0.377945227    5.669178402
H    -6.236096242   -4.157397495   -2.078698747
"#
            ),
            Err("Invalid unit")
        );

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
au
O     0.188972613    0.377945227    5.669178402
"#
            ),
            Err("The number of atoms does not agree with the number of coordinates found")
        );
    }
}
