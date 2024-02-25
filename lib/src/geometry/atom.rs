use crate::misc::elements::Element;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub el: Element,
    pub origin: [f64; 3],
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

    pub fn scale_coord(&mut self, factor: &f64) {
        self.origin.iter_mut().for_each(|i| *i *= factor);
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{} {:>14.9} {:>14.9} {:>14.9}",
            self.el, self.origin[0], self.origin[1], self.origin[2]
        )
    }
}

impl FromStr for Atom {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() != 4 {
            return Err("Expected format: 'element x y z'");
        }

        let element = Element::from_str(parts[0])?;
        let x = parts[1].parse().map_err(|_| "Invalid float for x")?;
        let y = parts[2].parse().map_err(|_| "Invalid float for y")?;
        let z = parts[3].parse().map_err(|_| "Invalid float for z")?;

        Ok(Atom::new(element, [x, y, z]))
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::atom::Atom;
    use crate::misc::elements::Element;

    #[test]
    fn serialize_atom() {
        assert_eq!(
            Atom::new(Element::Sn, [0.1, -0.1, 5.0]).to_string(),
            "Sn    0.100000000   -0.100000000    5.000000000\n"
        );
    }

    #[test]
    fn deserialize_atom() {
        use std::str::FromStr;

        // passing deserialisation
        assert_eq!(
            Atom::from_str("Sn    0.100000000   -0.100000000    5.000000000"),
            Ok(Atom::new(Element::Sn, [0.1, -0.1, 5.0]))
        );

        // failing deserialisations
        assert_eq!(
            Atom::from_str("Sn    0.100000000   -0.100000000"),
            Err("Expected format: 'element x y z'")
        );

        assert_eq!(
            Atom::from_str("Sn    o.100000000   -0.100000000    5.000000000"),
            Err("Invalid float for x")
        );

        assert_eq!(
            Atom::from_str("Sn    0.100000000   -0.o00000000    5.000000000"),
            Err("Invalid float for y")
        );

        assert_eq!(
            Atom::from_str("Sn    0.100000000   -0.100000000    5,000000000"),
            Err("Invalid float for z")
        );
    }
}
