use crate::tools::math::r#trait::MyTrait;
use std::ops::{Index, IndexMut, Mul};

pub type FVector = BaseVector<f64>;
pub type IVector = BaseVector<i64>;

#[derive(Clone, Default)]
pub struct BaseVector<T>
where
    T: MyTrait,
{
    n: usize,
    elem: Vec<T>,
}

// simple getter
impl<T: MyTrait> BaseVector<T> {
    pub fn size(&self) -> usize {
        self.n
    }
}

impl<T: MyTrait> BaseVector<T> {
    pub fn new(n: &usize) -> Self {
        if *n <= 0usize {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        BaseVector::<T> {
            n: *n,
            elem: Vec::with_capacity(*n),
        }
    }

    pub fn init(&mut self, value: &T) {
        self.elem = vec![value.clone(); self.n];
    }

    pub fn print(&self) {
        for i in 0..self.n {
            println!("{} ", self[i]);
        }
    }
}

impl<T: MyTrait> Index<usize> for BaseVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elem[index]
    }
}

impl<T: MyTrait> IndexMut<usize> for BaseVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elem[index]
    }
}

// allows FVector * scalar
impl Mul<f64> for FVector {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self::Output {
        Self {
            n: self.n,
            elem: self.elem.iter().map(|v| v * scalar).collect(),
        }
    }
}

// allows scalar * FVector
impl Mul<FVector> for f64 {
    type Output = FVector;

    fn mul(self, rhs: FVector) -> Self::Output {
        FVector {
            n: rhs.n,
            elem: rhs.elem.iter().map(|v| v * self).collect(),
        }
    }
}
