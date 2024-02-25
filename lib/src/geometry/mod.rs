use std::{fmt::Display, str::FromStr};

pub mod atom;
pub mod molecule;

#[derive(Clone, Debug, PartialEq, Default)]
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
