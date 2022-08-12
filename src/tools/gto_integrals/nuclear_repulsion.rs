use crate::tools::geometry::{distance, Atom};

pub fn nuclear_repulsion(atoms: &Vec<Atom>) -> f64 {
    let n = atoms.len();
    let mut vnn = 0.0;
    for i in 0..n {
        for j in 0..i {
            vnn += atoms[i].charge() as f64 * atoms[j].charge() as f64
                / distance(&atoms[i], &atoms[j]);
        }
    }
    vnn
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tools::constants;

    #[test]
    fn nuclear_repulsion() {
        let o = Atom::new(8.0, 15.999, [0.0000000000, 0.0000000000, -0.1190150726]);
        let h1 = Atom::new(1.0, 1.008, [0.7685504811, 0.0000000000, 0.4760602904]);
        let h2 = Atom::new(1.0, 1.008, [-0.7685504811, 0.0000000000, 0.4760602904]);

        let mut atoms = vec![o, h1, h2];
        for i in 0..atoms.len() {
            atoms[i].scale_coords(&constants::_ANG_AU);
        }

        assert_eq!(super::nuclear_repulsion(&atoms), 9.055003146181436);
    }
}
