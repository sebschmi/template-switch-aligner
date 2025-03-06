use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

use num_traits::{Bounded, CheckedAdd, CheckedSub, SaturatingSub, Zero};

/// The cost of an A* node.
pub trait AStarCost:
    From<Self::CostType>
    + From<u8>
    + Add<Output = Self>
    + Sub<Output = Self>
    + SaturatingSub
    + CheckedAdd
    + CheckedSub
    + AddAssign
    + SubAssign
    + FromStr
    + Zero
    + Bounded
    + Display
    + Debug
    + Ord
    + Eq
    + Hash
    + Copy
{
    type CostType;

    fn as_f64(&self) -> f64;

    fn as_u64(&self) -> u64;

    fn as_primitive(&self) -> Self::CostType;
}

macro_rules! primitive_cost {
    ($name:ident, $primitive:ident) => {
        #[doc = concat!("The cost of an A* node.\n\nThis cost type uses [`", stringify!($primitive), "`] as the internal representation of cost.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name($primitive);

        impl $crate::cost::AStarCost for $name {
            type CostType = $primitive;

            fn as_f64(&self) -> f64 {
                self.0 as f64
            }

            fn as_u64(&self)->u64{
                self.0 as u64
            }

            fn as_primitive(&self) -> Self::CostType {
                self.0
            }
        }

        impl From<$primitive> for $name {
            fn from(value: $primitive) -> Self {
                Self(value)
            }
        }

        impl From<u8> for $name {
            fn from(value: u8) -> Self {
                Self(value.into())
            }
        }

        impl std::ops::Add for $name {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                Self(self.0.checked_add(rhs.0).unwrap())
            }
        }

        impl std::ops::Sub for $name {
            type Output = Self;

            fn sub(self, rhs: Self) -> Self::Output {
                Self(self.0.checked_sub(rhs.0).unwrap())
            }
        }

        impl num_traits::SaturatingSub for $name {
            fn saturating_sub(&self, rhs: &Self) -> Self {
                Self(self.0.saturating_sub(rhs.0))
            }
        }

        impl num_traits::CheckedAdd for $name {
            fn checked_add(&self, rhs: &Self) -> Option<Self> {
                self.0.checked_add(rhs.0).map($name)
            }
        }

        impl num_traits::CheckedSub for $name {
            fn checked_sub(&self, rhs: &Self) -> Option<Self> {
                self.0.checked_sub(rhs.0).map($name)
            }
        }

        impl std::ops::AddAssign for $name {
            fn add_assign(&mut self, rhs: Self) {
                *self = *self + rhs;
            }
        }

        impl std::ops::SubAssign for $name {
            fn sub_assign(&mut self, rhs: Self) {
                *self = *self - rhs;
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                std::fmt::Display::fmt(&self.0, f)
            }
        }

        impl std::str::FromStr for $name {
            type Err = <$primitive as std::str::FromStr>::Err;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                $primitive::from_str(s).map(Self)
            }
        }

        impl num_traits::Zero for $name {
            fn zero() -> Self {
                $primitive::zero().into()
            }

            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }

        impl num_traits::Bounded for $name {
            fn min_value() -> Self {
                <$primitive as num_traits::Bounded>::min_value().into()
            }

            fn max_value() -> Self {
                <$primitive as num_traits::Bounded>::max_value().into()
            }
        }
    };
}

primitive_cost!(U16Cost, u16);
primitive_cost!(I16Cost, i16);
primitive_cost!(U32Cost, u32);
primitive_cost!(I32Cost, i32);
primitive_cost!(U64Cost, u64);
primitive_cost!(I64Cost, i64);
