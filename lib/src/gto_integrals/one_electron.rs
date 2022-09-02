use crate::basis::{cdim, dim, BasisShell, CartesianBasisFunction};
use crate::gto_integrals::integral_interface::IntegralInterface;
use crate::gto_integrals::kinetic_energy::kinetic_energy;
use crate::gto_integrals::nuclear_electron_attraction::nuclear_electron_attraction;
use crate::gto_integrals::overlap::overlap;
use crate::math::matrix::FMatrix;
use crate::math::matrix_symmetric::FMatrixSym;

pub enum OneElectronKernel {
    Kinetic,
    NuclearAttraction,
    Overlap,
}

impl IntegralInterface {
    fn cartesian_to_spherical_transformation(
        &self,
        l_a: &u8,
        l_b: &u8,
        matrix_cartesian: FMatrix,
    ) -> FMatrix {
        let dim_a = dim(l_a);
        let dim_b = dim(l_b);
        let cdim_a = cdim(l_a);
        let cdim_b = cdim(l_b);

        let ta = self.basis().trafo_matrix(l_a);
        let tb = self.basis().trafo_matrix(l_b);
        let mut matrix_spherical = FMatrix::zero(dim_a, dim_b);

        for i in 0..dim_a {
            for j in 0..dim_b {
                for a in 0..cdim_a {
                    for b in 0..cdim_b {
                        matrix_spherical[(i, j)] +=
                            ta[(i, a)] * tb[(j, b)] * matrix_cartesian[(a, b)]
                    }
                }
            }
        }

        matrix_spherical
    }

    fn calc_one_electron_cbf_cbf(
        &self,
        kernel: &OneElectronKernel,
        a: &CartesianBasisFunction,
        b: &CartesianBasisFunction,
    ) -> f64 {
        let mut value = 0.0;
        for (ia, ca) in a.coefs().iter().enumerate() {
            for (ib, cb) in b.coefs().iter().enumerate() {
                let integral = {
                    match kernel {
                        OneElectronKernel::Kinetic => kinetic_energy(
                            &a.exps()[ia],
                            &a.ml_i16(),
                            a.origin(),
                            &b.exps()[ib],
                            &b.ml_i16(),
                            b.origin(),
                        ),
                        OneElectronKernel::NuclearAttraction => nuclear_electron_attraction(
                            &a.exps()[ia],
                            &a.ml_i16(),
                            a.origin(),
                            &b.exps()[ib],
                            &b.ml_i16(),
                            b.origin(),
                            self.atoms(),
                        ),
                        OneElectronKernel::Overlap => overlap(
                            &a.exps()[ia],
                            &a.ml_i16(),
                            a.origin(),
                            &b.exps()[ib],
                            &b.ml_i16(),
                            b.origin(),
                        ),
                    }
                };

                value += a.norm()[ia] * b.norm()[ib] * ca * cb * integral;
            }
        }

        value
    }

    fn calc_one_electron_shell_shell(
        &self,
        kernel: &OneElectronKernel,
        a: &BasisShell,
        b: &BasisShell,
    ) -> FMatrix {
        let dim_a = a.cdim();
        let dim_b = b.cdim();

        let mut matrix_cartesian = FMatrix::zero(dim_a, dim_b);

        for i in 0..dim_a {
            for j in 0..dim_b {
                matrix_cartesian[(i, j)] =
                    self.calc_one_electron_cbf_cbf(&kernel, &a.cbf()[i], &b.cbf()[j]);
            }
        }

        self.cartesian_to_spherical_transformation(a.l(), b.l(), matrix_cartesian)
    }

    pub fn calc_one_electron_integral(&self, kernel: OneElectronKernel) -> FMatrixSym {
        let dim = self.basis().dim();
        let mut one_electron_integral = FMatrixSym::zero(dim);

        for i in 0..self.basis().shells().len() {
            for j in 0..=i {
                let one_electron_sub_matrix = self.calc_one_electron_shell_shell(
                    &kernel,
                    &self.basis().shells()[i],
                    &self.basis().shells()[j],
                );
                // println!("{}", one_electron_sub_matrix);
                let offset_i = self.basis().offset(i);
                let offset_j = self.basis().offset(j);

                // todo: put into symmetric matrix (make function)
                for a in 0..one_electron_sub_matrix.rows() {
                    for b in 0..one_electron_sub_matrix.cols() {
                        one_electron_integral[(a + offset_i, b + offset_j)] =
                            one_electron_sub_matrix[(a, b)];
                    }
                }
            }
        }

        one_electron_integral
    }
}
