pub trait Factorial {
    /// The resulting type of the factorial functions
    type Output;

    /// Calculates the factorial of a number (n!)
    ///
    /// # Example
    /// ```
    /// // assert_eq!(5_u8.factorial(), 120);
    /// ```
    fn factorial(&self) -> Self::Output;

    /// Calculates the double factorial of a number (n!!)
    ///
    /// # Example
    /// ```
    /// // assert_eq!(5_u8.factorial2(), 15);
    /// ```
    fn factorial2(&self) -> Self::Output;
}

/// Type specific implementation for the Factorial trait using a macro
macro_rules! fact_unsigned_impl {
    ($($t:ty),*) => {$(

impl Factorial for $t {
    type Output = usize;

    fn factorial(&self) -> Self::Output {
        let n = *self as Self::Output;
        match self {
            0 | 1 => 1,
            _ => (1..=n).product(),
        }
    }

    fn factorial2(&self) -> Self::Output {
        let n = *self as Self::Output;
        match self {
            0 | 1 => 1,
            _ => (1..=n).rev().step_by(2).product()
        }
    }
}

    )*};
}

fact_unsigned_impl!(u8, u16, u32, u64, usize);

/// Type specific implementation for the Factorial trait using a macro
macro_rules! fact_signed_impl {
    ($($t:ty),*) => {$(

impl Factorial for $t {
    type Output = usize;

    fn factorial(&self) -> Self::Output {
        if (self <= &1) {
            return 1;
        }

        let n = *self as Self::Output;
        (1..=n).product()
    }

    fn factorial2(&self) -> Self::Output {
        if (self <= &1) {
            return 1;
        }

        let n = *self as Self::Output;
        (1..=n).rev().step_by(2).product()
    }
}

    )*};
}

fact_signed_impl!(i8, i16, i32, i64, isize);

pub trait BinomialCoefficient {
    type Output;

    /// Calculates the binomial coefficient
    /// as $(3 / 2)$
    fn binomial_coefficient(&self, k: &Self) -> Self::Output;
}

/// Type specific implementation for the Factorial trait using a macro
macro_rules! binomial_impl {
    ($($t:ty),*) => {$(

impl BinomialCoefficient for $t {
    type Output = usize;

    fn binomial_coefficient(&self, k: &Self) -> Self::Output {
        if k > self {
            panic!("k is larger than n!");
        }

        (1..=*k).map(|i| (*self + 1 - i) / i).product::<$t>() as usize
    }
}

    )*};
}

binomial_impl!(u8, u16, u32, u64, usize);

// todo: add test
pub fn gaussian_product_center(
    a: &f64,
    a_center: &[f64; 3],
    b: &f64,
    b_center: &[f64; 3],
) -> [f64; 3] {
    let mut product = [0.0; 3];
    let f = a + b;
    product[0] = (a * a_center[0] + b * b_center[0]) / f;
    product[1] = (a * a_center[1] + b * b_center[1]) / f;
    product[2] = (a * a_center[2] + b * b_center[2]) / f;

    product
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial() {
        // todo: add large cases
        assert_eq!(1, 0_u8.factorial());
        assert_eq!(1, 1_u8.factorial());
        assert_eq!(2, 2_u8.factorial());
        assert_eq!(6, 3_u8.factorial());
        assert_eq!(24, 4_u8.factorial());
        assert_eq!(120, 5_u8.factorial());
    }

    #[test]
    fn test_factorial2() {
        assert_eq!(1, 0_u8.factorial2());
        assert_eq!(1, 1_u8.factorial2());
        assert_eq!(2, 2_u8.factorial2());
        assert_eq!(3, 3_u8.factorial2());
        assert_eq!(8, 4_u8.factorial2());
        assert_eq!(15, 5_u8.factorial2());
        assert_eq!(48, 6_u8.factorial2());
    }

    #[test]
    fn test_binomial_coeff() {
        assert_eq!(10, 5_u8.binomial_coefficient(&2));
    }
}
