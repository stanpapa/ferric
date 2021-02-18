use rust_scf::tools::*;
// use crate::basis::CartesianBasisFunction;

fn main() {

  let o = geometry::Atom  { charge: 8.0, mass: 15.999, origin: [ 0.0000000000, 0.0000000000,-0.1190150726]};
  let h1 = geometry::Atom { charge: 1.0, mass: 1.008,  origin: [ 0.7685504811, 0.0000000000, 0.4760602904]};
  let h2 = geometry::Atom { charge: 1.0, mass: 1.008,  origin: [-0.7685504811, 0.0000000000, 0.4760602904]};

  let mut atoms = vec![o, h1, h2];
  for i in 0..atoms.len() {
    atoms[i].scale_coords(&crate::constants::_ANG_AU);
  }
  println!("Vnn:      {}", crate::gto_integrals::nuclear_repulsion::nuclear_repulsion(&atoms));

  crate::gto_bases::def2_tzvp::load_def2_tzvp();
  
  let imat = crate::math::matrix::IMatrix::new(&2);
  // println!("{:?}", imat);

}
