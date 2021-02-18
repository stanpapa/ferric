use crate::tools::geometry::{Atom, distance};

pub fn nuclear_repulsion(atoms: &Vec<Atom>) -> f64 {
  
  let n = atoms.len();
  let mut vnn = 0.0;
  for i in 0..n {
    for j in 0..i {
      vnn += atoms[i].charge * atoms[j].charge / distance(&atoms[i], &atoms[j]);
    }
  }
  vnn
}
