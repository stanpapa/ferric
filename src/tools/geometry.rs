pub struct Atom {
  pub charge: f64,
  pub mass: f64,
  pub origin: [f64; 3],
}

impl Atom {

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

  (dx*dx + dy*dy + dz*dz).sqrt()
}
