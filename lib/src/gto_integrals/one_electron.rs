use crate::{
    gto_basis_sets::basis::{BasisShell, CartesianBasisFunction},
    gto_integrals::{
        h_core::h_core, integral_interface::IntegralInterface, kinetic_energy::kinetic_energy,
        nuclear_electron_attraction::nuclear_electron_attraction, overlap::overlap,
    },
    linear_algebra::matrix::FMatrix,
};

use std::{fmt::Display, slice::Iter, time::Instant};

pub enum OneElectronKernel {
    HCore,
    Kinetic,
    NuclearAttraction,
    Overlap,
}

impl OneElectronKernel {
    pub fn to_filename(&self) -> &str {
        match self {
            OneElectronKernel::HCore => "h_ao.tmp",
            OneElectronKernel::Kinetic => "t_ao.tmp",
            OneElectronKernel::NuclearAttraction => "v_ao.tmp",
            OneElectronKernel::Overlap => "s_ao.tmp",
        }
    }

    pub fn iter() -> Iter<'static, OneElectronKernel> {
        static KERNEL: [OneElectronKernel; 4] = [
            OneElectronKernel::HCore,
            OneElectronKernel::Kinetic,
            OneElectronKernel::NuclearAttraction,
            OneElectronKernel::Overlap,
        ];
        KERNEL.iter()
    }
}

impl Display for OneElectronKernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OneElectronKernel::HCore => write!(f, "H")?,
            OneElectronKernel::Kinetic => write!(f, "T")?,
            OneElectronKernel::NuclearAttraction => write!(f, "V")?,
            OneElectronKernel::Overlap => write!(f, "S")?,
        }

        Ok(())
    }
}

impl IntegralInterface {
    fn cartesian_to_spherical_transformation_1e(
        &self,
        l_a: &u8,
        l_b: &u8,
        matrix_cartesian: FMatrix,
    ) -> FMatrix {
        let ta = self.basis().trafo_matrix(l_a);
        let tb = self.basis().trafo_matrix(l_b);
        ta * (matrix_cartesian * tb.transposed())
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
                        OneElectronKernel::HCore => h_core(
                            &a.exps()[ia],
                            &a.ml_i16(),
                            a.origin(),
                            &b.exps()[ib],
                            &b.ml_i16(),
                            b.origin(),
                            self.atoms(),
                        ),
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
                    self.calc_one_electron_cbf_cbf(kernel, &a.cbf()[i], &b.cbf()[j]);
            }
        }

        self.cartesian_to_spherical_transformation_1e(a.l(), b.l(), matrix_cartesian)
    }

    pub fn calc_one_electron_integral(&self, kernel: OneElectronKernel) -> FMatrix {
        print!("Calculating one-electron integral: {} ... ", kernel);
        let t = Instant::now();

        let dim = self.basis().dim();
        let mut one_electron_integral = FMatrix::zero(dim, dim);

        for i in 0..self.basis().shells().len() {
            for j in 0..=i {
                let one_electron_sub_matrix = self.calc_one_electron_shell_shell(
                    &kernel,
                    &self.basis().shells()[i],
                    &self.basis().shells()[j],
                );
                let offset_i = self.basis().offset(i);
                let offset_j = self.basis().offset(j);

                for a in 0..one_electron_sub_matrix.rows {
                    for b in 0..one_electron_sub_matrix.cols {
                        one_electron_integral[(a + offset_i, b + offset_j)] =
                            one_electron_sub_matrix[(a, b)];
                    }
                }
            }
        }

        // fill opposite triangle
        for i in 0..one_electron_integral.rows {
            for j in 0..i {
                one_electron_integral[(j, i)] = one_electron_integral[(i, j)];
            }
        }

        println!("done ({:?})", t.elapsed());

        one_electron_integral.store(kernel.to_filename());
        one_electron_integral
    }
}
