use crate::misc::elements::Element;
use std::fmt::{Display, Formatter};
use std::str::FromStr;

use serde::{Deserialize, Serialize};
// use toml::ser::Serializer; // todo! Derive (De)Serialize myself?

use toml::value::Value;

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq)]
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
            "{} {:>14.9} {:>14.9} {:>14.9}",
            self.el, self.origin[0], self.origin[1], self.origin[2]
        )
    }
}

impl Atom {
    pub fn from_array(value: &Value) -> Self {
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

// impl Serialize for Atom {
//     fn serialize<S>(&self, serializer: S) -> Result<serde::ser::Ok, serde::ser::Error>
//     where
//         S: Serializer,
//     {
//         let toml
//         Ok(serde::ser::Ok)
//     }
// }

#[cfg(test)]
mod tests {
    use crate::geometry::atom::Atom;
    use crate::misc::elements::Element;

    #[test]
    fn serialize_atom() {
        assert_eq!(
            toml::to_string(&Atom::new(Element::Bi, [0.1, -0.1, 5.0])).unwrap(),
            r#"el = "Bi"
origin = [0.1, -0.1, 5.0]
"#
        );
    }

    #[test]
    fn deserialize_atom() {
        assert_eq!(
            toml::from_str::<Atom>(
                r#"el = "Bi"
origin = [0.1, -0.1, 5.0]
"#
            )
            .unwrap(),
            Atom::new(Element::Bi, [0.1, -0.1, 5.0])
        );
    }
}
