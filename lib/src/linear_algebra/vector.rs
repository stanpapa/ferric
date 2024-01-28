use crate::linear_algebra::scalar::Scalar;
use std::fmt::{Display, Formatter};
use std::ops::{Deref, DerefMut, Index, IndexMut};

pub type FVector = Vector<f64>;
pub type IVector = Vector<i64>;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Vector<T: Scalar> {
    pub n: usize,
    data: Vec<T>,
}

// simple getter
impl<T: Scalar> Vector<T> {
    fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_ok(&self) -> bool {
        self.n == self.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<T: Scalar> Vector<T> {
    pub fn new(n: usize) -> Self {
        if n <= 0 {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        Vector::<T> {
            n,
            data: Vec::with_capacity(n),
        }
    }

    pub fn new_with_value(n: usize, value: T) -> Self {
        if n <= 0 {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        Vector::<T> {
            n,
            data: vec![value; n],
        }
    }

    pub fn new_from_vec(v: &[T]) -> Self {
        Self {
            n: v.len(),
            data: v.to_vec(),
        }
    }

    pub fn zero(n: usize) -> Self {
        Self::new_with_value(n, T::default())
    }

    pub fn init(&mut self, value: T) {
        self.data = vec![value; self.n];
    }

    // fn del(&mut self) {
    //     self.n = 0;
    //     self.elem.clear();
    // }

    pub fn print(&self) {
        if self.is_empty() {
            println!("Empty vector. Nothing to print.");
            return;
        }

        for i in 0..self.n {
            println!("{} ", self[i]);
        }
    }
}

// implement `Deref` so, Vector<T>.iter() can be called
impl<T: Scalar> Deref for Vector<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

// implement `DerefMut` so, Vector<T>.iter_mut() can be called
impl<T: Scalar> DerefMut for Vector<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

// ---------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------

impl<T: Scalar> Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T: Scalar> IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<T: Scalar> Display for Vector<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for i in &self.data {
            writeln!(f, "{:14.9}", i)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 5;

    #[test]
    fn new() {
        let vec = FVector {
            n: N,
            data: Vec::with_capacity(N),
        };

        assert_eq!(vec, FVector::new(N));
    }

    #[test]
    fn new_with_value() {
        let vec = FVector {
            n: N,
            data: vec![5.0; N],
        };

        assert_eq!(vec, FVector::new_with_value(N, 5.0));
    }

    #[test]
    fn new_from_vec() {
        let vec = IVector {
            n: N,
            data: vec![0, 1, 2, 3, 4],
        };

        assert_eq!(vec, IVector::new_from_vec(&[0, 1, 2, 3, 4]));
    }

    #[test]
    fn zero() {
        let vec = FVector {
            n: N,
            data: vec![0.0; N],
        };

        assert_eq!(vec, FVector::zero(N));
    }

    #[test]
    fn init() {
        let mut vec = FVector::new(N);
        vec.init(5.0);

        assert_eq!(vec, FVector::new_with_value(N, 5.0));
    }

    #[test]
    fn index() {
        let vec = IVector {
            n: N,
            data: (0..5).collect(),
        };

        assert_eq!(vec[3], 3);
    }
}
