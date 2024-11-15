use std::{
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

use num_traits::{CheckedSub, SaturatingSub};

type CostType = u64;

/// The cost of an alignment.
///
/// This cost type is not allowed to be negative.
/// This is important for example when using Dijkstra or A* to compute an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cost(CostType);

impl Cost {
    pub const MIN: Self = Self(CostType::MIN);
    pub const MAX: Self = Self(CostType::MAX);
    pub const ZERO: Self = Self(0);

    pub fn as_u64(&self) -> u64 {
        self.0
    }
}

impl From<CostType> for Cost {
    fn from(value: CostType) -> Self {
        Self(value)
    }
}

impl Add for Cost {
    type Output = Cost;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Cost {
    type Output = Cost;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl SaturatingSub for Cost {
    fn saturating_sub(&self, rhs: &Self) -> Self {
        Self(self.0.saturating_sub(rhs.0))
    }
}

impl CheckedSub for Cost {
    fn checked_sub(&self, rhs: &Self) -> Option<Self> {
        self.0.checked_sub(rhs.0).map(Cost)
    }
}

impl AddAssign for Cost {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Cost {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl std::fmt::Display for Cost {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl FromStr for Cost {
    type Err = <CostType as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        CostType::from_str(s).map(Self)
    }
}
