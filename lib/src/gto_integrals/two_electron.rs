use crate::{
    gto_basis_sets::basis::{cdim, dim, BasisShell, CartesianBasisFunction},
    gto_integrals::{eri::eri, integral_interface::IntegralInterface},
    linear_algebra::{
        matrix::FMatrix, matrix_container::FMatrixContainer, matrix_container::FMatrixSymContainer,
        matrix_symmetric::FMatrixSym,
    },
};

use std::slice::Iter;

pub enum TwoElectronKernel {
    ERI,
}

impl TwoElectronKernel {
    pub fn to_filename(&self) -> &str {
        match self {
            TwoElectronKernel::ERI => "eri_ao.tmp",
        }
    }

    pub fn iter() -> Iter<'static, TwoElectronKernel> {
        static KERNEL: [TwoElectronKernel; 1] = [TwoElectronKernel::ERI];
        KERNEL.iter()
    }
}

impl IntegralInterface {
    // pub fn calc_two_electron_integral(&self, kernel: TwoElectronKernel) {
    pub fn calc_two_electron_integral(&self, kernel: TwoElectronKernel) -> FMatrixSymContainer {
        let dim = self.basis().dim();
        let mut two_electron_integral = FMatrixSymContainer::new();

        for i in 0..self.basis().shells().len() {
            for j in 0..=i {
                for k in 0..self.basis().shells().len() {
                    for l in 0..=k {
                        let integral_sub = self.calc_two_electron_shell(
                            &kernel,
                            &self.basis().shells()[i],
                            &self.basis().shells()[j],
                            &self.basis().shells()[k],
                            &self.basis().shells()[l],
                        );

                        let offset_i = self.basis().offset(i);
                        let offset_j = self.basis().offset(j);
                        let offset_k = self.basis().offset(k);
                        let offset_l = self.basis().offset(l);

                        for ab in integral_sub.keys() {
                            let ab_offset = (ab.0 + offset_i, ab.1 + offset_j);

                            // only keep lower triangle
                            if ab_offset.0 < ab_offset.1 {
                                continue;
                            }

                            let tmp = two_electron_integral
                                .entry(&ab_offset)
                                .or_insert_with(|| FMatrixSym::zero(dim));

                            // an entry exists at ab, so it can be safely retrieved
                            let sub = &integral_sub[*ab];
                            for c in 0..sub.rows {
                                for d in 0..sub.cols {
                                    tmp[(c + offset_k, d + offset_l)] = sub[(c, d)];
                                }
                            }
                        }
                    }
                }
            }
        }

        // two_electron_integral.store(kernel.to_filename());
        two_electron_integral
    }

    fn calc_two_electron_shell(
        &self,
        kernel: &TwoElectronKernel,
        a: &BasisShell,
        b: &BasisShell,
        c: &BasisShell,
        d: &BasisShell,
    ) -> FMatrixContainer {
        let dim_a = a.cdim();
        let dim_b = b.cdim();
        let dim_c = c.cdim();
        let dim_d = d.cdim();

        let mut integral_cartesian = FMatrixContainer::new();
        for i in 0..dim_a {
            for j in 0..dim_b {
                let mut mat = FMatrix::zero(dim_c, dim_d);
                for k in 0..dim_c {
                    for l in 0..dim_d {
                        mat[(k, l)] = self.calc_two_electron_cbf(
                            kernel,
                            &a.cbf()[i],
                            &b.cbf()[j],
                            &c.cbf()[k],
                            &d.cbf()[l],
                        );
                    }
                }
                integral_cartesian.insert((i, j), &mat);
            }
        }

        // transform to spherical
        self.cartesian_to_spherical_transformation_2e(
            a.l(),
            b.l(),
            c.l(),
            d.l(),
            integral_cartesian,
        )
    }

    fn calc_two_electron_cbf(
        &self,
        kernel: &TwoElectronKernel,
        a: &CartesianBasisFunction,
        b: &CartesianBasisFunction,
        c: &CartesianBasisFunction,
        d: &CartesianBasisFunction,
    ) -> f64 {
        let mut value = 0.0;
        for (ia, ca) in a.coefs().iter().enumerate() {
            for (ib, cb) in b.coefs().iter().enumerate() {
                for (ic, cc) in c.coefs().iter().enumerate() {
                    for (id, cd) in d.coefs().iter().enumerate() {
                        let integral = match kernel {
                            TwoElectronKernel::ERI => eri(
                                &a.exps()[ia],
                                &a.ml_i16(),
                                a.origin(),
                                &b.exps()[ib],
                                &b.ml_i16(),
                                b.origin(),
                                &c.exps()[ic],
                                &c.ml_i16(),
                                c.origin(),
                                &d.exps()[id],
                                &d.ml_i16(),
                                d.origin(),
                            ),
                        };

                        value += a.norm()[ia]
                            * b.norm()[ib]
                            * c.norm()[ic]
                            * d.norm()[id]
                            * ca
                            * cb
                            * cc
                            * cd
                            * integral;
                    }
                }
            }
        }

        value
    }

    fn cartesian_to_spherical_transformation_2e(
        &self,
        l_a: &u8,
        l_b: &u8,
        l_c: &u8,
        l_d: &u8,
        cartesian: FMatrixContainer,
    ) -> FMatrixContainer {
        let ta = self.basis().trafo_matrix(l_a);
        let tb = self.basis().trafo_matrix(l_b);
        let tc = self.basis().trafo_matrix(l_c);
        let td = self.basis().trafo_matrix(l_d);

        let mut spherical = FMatrixContainer::new();
        // todo: use BLAS
        for i in 0..dim(l_a) {
            for j in 0..dim(l_b) {
                spherical.insert((i, j), &FMatrix::zero(dim(l_c), dim(l_d)));
                for k in 0..dim(l_c) {
                    for l in 0..dim(l_d) {
                        for a in 0..cdim(l_a) {
                            for b in 0..cdim(l_b) {
                                for c in 0..cdim(l_c) {
                                    for d in 0..cdim(l_d) {
                                        spherical[(i, j)][(k, l)] += ta[(i, a)]
                                            * tb[(j, b)]
                                            * tc[(k, c)]
                                            * td[(l, d)]
                                            * cartesian[(a, b)][(c, d)];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        spherical
    }
}
