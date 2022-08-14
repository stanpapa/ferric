use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

pub trait Scalar:
    Copy
    + Clone
    + Display
    + Default
    + Add
    + AddAssign
    + Div
    + DivAssign
    + Mul
    + MulAssign
    + Sub
    + SubAssign
    + PartialEq
{
}

impl<
        T: Copy
            + Clone
            + Display
            + Default
            + Add
            + AddAssign
            + Div
            + DivAssign
            + Mul
            + MulAssign
            + Sub
            + SubAssign
            + PartialEq,
    > Scalar for T
{
}

// /// The basic scalar type for all structures of `...`.
// ///
// /// This does not make any assumption on the algebraic properties of `Self`.
// /// Taken from `nalgebra`
// pub trait ScalarNew: 'static + Clone + PartialEq + Debug {}
//
// impl<T: 'static + Clone + PartialEq + Debug> ScalarNew for T {}
