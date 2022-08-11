use crate::tools::math::scalar::Scalar;
use std::iter::FromIterator;
use std::ops::{Index, IndexMut, Mul, MulAssign};

pub type FVector = BaseVector<f64>;
pub type IVector = BaseVector<i64>;

#[derive(Clone, Default)]
pub struct BaseVector<T: Scalar> {
    n: usize,
    elem: Vec<T>,
}

// simple getter
impl<T: Scalar> BaseVector<T> {
    pub fn size(&self) -> usize {
        self.n
    }
}

impl<T: Scalar> BaseVector<T> {
    pub fn new(n: usize) -> Self {
        if n <= 0 {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        BaseVector::<T> {
            n,
            elem: vec![T::default(); n],
        }
    }

    pub fn init(&mut self, value: T) {
        self.elem = vec![value; self.n];
    }

    pub fn print(&self) {
        for i in 0..self.n {
            println!("{} ", self[i]);
        }
    }
}

impl<T: Scalar> Index<usize> for BaseVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elem[index]
    }
}

impl<T: Scalar> IndexMut<usize> for BaseVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elem[index]
    }
}

// allows BaseVector * scalar
impl<T: Scalar> Mul<T> for BaseVector<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseVector<T>;

    fn mul(self, scalar: T) -> Self::Output {
        Self {
            elem: self.elem.into_iter().map(|v| v * scalar).collect(),
            ..self
        }
    }
}

// allows scalar * BaseVector
// todo: make it work with generics
// impl<T: DataTrait> Mul<BaseVector<T>> for T {
//     type Output = BaseVector<T>;
//
//     fn mul(self, rhs: BaseVector<T>) -> Self::Output {
//         BaseVector::<T> {
//             n: rhs.n,
//             elem: rhs.elem.into_iter().map(|v| v * self).collect(),
//         }
//     }
// }

// allows BaseMatrix *= scalar
impl<T: Scalar> MulAssign<T> for BaseVector<T> {
    fn mul_assign(&mut self, rhs: T) {
        for x in &mut self.elem {
            *x *= rhs;
        }
    }
}
// todo: add tests
