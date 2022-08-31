use crate::elements::Element;
use crate::elements::Element::{H, O};

pub struct Molecule {
    atoms: Vec<Atom>,
    charge: i8,
    multiplicity: u8,
    num_electrons_alpha: usize,
    num_electrons_beta: usize,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>, charge: i8, multiplicity: u8) -> Molecule {
        let mut num_electrons: i8 = -charge;
        num_electrons += atoms.iter().map(|a| a.z()).sum::<u8>() as i8;

        if num_electrons <= 0 {
            panic!("Invalid charge of {}.", charge);
        }

        // Mult = 2Ms + 1 thus the number of unpaired electrons (taken as α) is Mult-1 = 2Ms
        let alpha_excess = multiplicity as usize - 1;

        // the number of β electrons must be equal the number of doubly occupied orbitals
        // Ndo = (nelec - n_of_unpaired_elec)/2 this must be an integer
        // if isodd(num_electrons) != isodd(alpha_excess) {
        if (num_electrons % 2 == 0) && (alpha_excess % 2 == 1) {
            panic!(
                "Incompatible charge {} and multiplicity {}",
                charge, multiplicity
            );
        }

        let num_electrons_beta = (num_electrons as usize - alpha_excess) / 2;
        let num_electrons_alpha = num_electrons as usize - num_electrons_beta;

        Molecule {
            atoms,
            charge,
            multiplicity,
            num_electrons_alpha,
            num_electrons_beta,
        }
    }

    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &Vec<Atom> {
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
        Self::new(
            vec![
                Atom::new(O, [0.0000000000, 0.0000000000, -0.1190150726]),
                Atom::new(H, [0.7685504811, 0.0000000000, 0.4760602904]),
                Atom::new(H, [-0.7685504811, 0.0000000000, 0.4760602904]),
            ],
            0,
            1,
        )
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
        for i in 0..3 {
            self.origin[i] *= factor;
        }
    }
}

pub fn distance(a: &Atom, b: &Atom) -> f64 {
    let dx = a.origin[0] - b.origin[0];
    let dy = a.origin[1] - b.origin[1];
    let dz = a.origin[2] - b.origin[2];

    (dx * dx + dy * dy + dz * dz).sqrt()
}
