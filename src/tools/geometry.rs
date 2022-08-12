#[derive(Copy, Clone)]
pub struct Atom {
    charge: f32,
    mass: f32,
    origin: [f64; 3],
}

impl Atom {
    pub fn new(charge: f32, mass: f32, origin: [f64; 3]) -> Self {
        Self {
            charge,
            mass,
            origin,
        }
    }

    pub fn scale_coords(&mut self, factor: &f64) {
        for i in 0..3 {
            self.origin[i] *= factor;
        }
    }

    pub fn charge(&self) -> f32 {
        self.charge
    }

    pub fn mass(&self) -> f32 {
        self.mass
    }

    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }
}

pub fn distance(a: &Atom, b: &Atom) -> f64 {
    let dx = a.origin[0] - b.origin[0];
    let dy = a.origin[1] - b.origin[1];
    let dz = a.origin[2] - b.origin[2];

    (dx * dx + dy * dy + dz * dz).sqrt()
}
