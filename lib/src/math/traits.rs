use crate::math::scalar::Scalar;
use crate::math::vector::Vector;

pub trait Dot<T: Scalar> {
    type Output;

    fn dot(&self, rhs: &Vector<T>) -> Self::Output;
}

pub trait Norm {
    type Output;

    fn norm(&self) -> Self::Output;
}
