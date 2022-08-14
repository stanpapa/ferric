use std::any::{Any, TypeId};
use std::fmt::Debug;
use std::ops::{Add, Sub};

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
struct Dynamic {
    value: usize,
}

impl Dynamic {
    /// A dynamic size equal to `value`.
    #[inline]
    pub const fn new(value: usize) -> Self {
        Self { value }
    }
}

/// Trait implemented by `Dynamic`.
pub trait IsDynamic {}
/// Trait implemented by `Dynamic` and type-level integers different from `U1`.
pub trait IsNotStaticOne {}

impl IsDynamic for Dynamic {}
impl IsNotStaticOne for Dynamic {}

/// Dim of dynamically-sized algebraic entities.
/// Trait implemented by any type that can be used as a dimension. This includes type-level
/// integers and `Dynamic` (for dimensions not known at compile-time).
pub trait Dim: Any + Debug + Copy + PartialEq {
    #[inline(always)]
    fn is<D: Dim>() -> bool {
        TypeId::of::<Self>() == TypeId::of::<D>()
    }

    /// Gets the compile-time value of `Self`. Returns `None` if it is not known, i.e., if `Self =
    /// Dynamic`.
    fn try_to_usize() -> Option<usize>;

    /// Gets the run-time value of `self`. For type-level integers, this is the same as
    /// `Self::try_to_usize().unwrap()`.
    fn value(&self) -> usize;

    /// Builds an instance of `Self` from a run-time value. Panics if `Self` is a type-level
    /// integer and `dim != Self::try_to_usize().unwrap()`.
    fn from_usize(dim: usize) -> Self;
}

impl Dim for Dynamic {
    #[inline]
    fn try_to_usize() -> Option<usize> {
        None
    }

    #[inline]
    fn value(&self) -> usize {
        self.value
    }

    #[inline]
    fn from_usize(dim: usize) -> Self {
        Self::new(dim)
    }
}

impl Add<usize> for Dynamic {
    type Output = Dynamic;

    #[inline]
    fn add(self, rhs: usize) -> Self {
        Self::new(self.value + rhs)
    }
}

impl Sub<usize> for Dynamic {
    type Output = Dynamic;

    #[inline]
    fn sub(self, rhs: usize) -> Self {
        Self::new(self.value - rhs)
    }
}

/*
 *
 * Operations.
 *
 */

macro_rules! dim_ops(
    ($($DimOp:    ident, $DimNameOp: ident,
       $Op:       ident, $op: ident, $op_path: path,
       $DimResOp: ident, $DimNameResOp: ident,
       $ResOp: ident);* $(;)*) => {$(
        pub type $DimResOp<D1, D2> = <D1 as $DimOp<D2>>::Output;

        pub trait $DimOp<D: Dim>: Dim {
            type Output: Dim;

            fn $op(self, other: D) -> Self::Output;
        }

        impl<const A: usize, const B: usize> $DimOp<Const<B>> for Const<A>
        where
            Const<A>: ToTypenum,
            Const<B>: ToTypenum,
            <Const<A> as ToTypenum>::Typenum: $Op<<Const<B> as ToTypenum>::Typenum>,
            $ResOp<<Const<A> as ToTypenum>::Typenum, <Const<B> as ToTypenum>::Typenum>: ToConst,
        {
            type Output =
                <$ResOp<<Const<A> as ToTypenum>::Typenum, <Const<B> as ToTypenum>::Typenum> as ToConst>::Const;

            fn $op(self, _: Const<B>) -> Self::Output {
                Self::Output::name()
            }
        }

        impl<D: Dim> $DimOp<D> for Dynamic {
            type Output = Dynamic;

            #[inline]
            fn $op(self, other: D) -> Dynamic {
                Dynamic::new($op_path(self.value, other.value()))
            }
        }

        // TODO: use Const<T> instead of D: DimName?
        impl<D: DimName> $DimOp<Dynamic> for D {
            type Output = Dynamic;

            #[inline]
            fn $op(self, other: Dynamic) -> Dynamic {
                Dynamic::new($op_path(self.value(), other.value))
            }
        }

        pub type $DimNameResOp<D1, D2> = <D1 as $DimNameOp<D2>>::Output;

        pub trait $DimNameOp<D: DimName>: DimName {
            type Output: DimName;

            fn $op(self, other: D) -> Self::Output;
        }

        impl<const A: usize, const B: usize> $DimNameOp<Const<B>> for Const<A>
        where
            Const<A>: ToTypenum,
            Const<B>: ToTypenum,
            <Const<A> as ToTypenum>::Typenum: $Op<<Const<B> as ToTypenum>::Typenum>,
            $ResOp<<Const<A> as ToTypenum>::Typenum, <Const<B> as ToTypenum>::Typenum>: ToConst,
        {
            type Output =
                <$ResOp<<Const<A> as ToTypenum>::Typenum, <Const<B> as ToTypenum>::Typenum> as ToConst>::Const;

            fn $op(self, _: Const<B>) -> Self::Output {
                Self::Output::name()
            }
        }
   )*}
);

dim_ops!(
    DimAdd, DimNameAdd, Add, add, Add::add, DimSum,     DimNameSum,     Sum;
    DimMul, DimNameMul, Mul, mul, Mul::mul, DimProd,    DimNameProd,    Prod;
    DimSub, DimNameSub, Sub, sub, Sub::sub, DimDiff,    DimNameDiff,    Diff;
    // DimDiv, DimNameDiv, Div, div, Div::div, DimQuot,    DimNameQuot,    Quot;
    DimMin, DimNameMin, Min, min, cmp::min, DimMinimum, DimNameMinimum, Minimum;
    DimMax, DimNameMax, Max, max, cmp::max, DimMaximum, DimNameMaximum, Maximum;
);
