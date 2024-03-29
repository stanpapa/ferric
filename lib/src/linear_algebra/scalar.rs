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
