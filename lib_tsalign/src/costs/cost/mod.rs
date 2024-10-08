use std::ops::{Add, AddAssign, Sub, SubAssign};

/// The cost of an alignment.
///
/// This cost type is not allowed to be negative.
/// This is important for example when using Dijkstra or A* to compute an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cost(u64);

impl Cost {
    pub const MAX: Self = Self(u64::MAX);
    pub const ZERO: Self = Self(0);

    pub fn as_u64(&self) -> u64 {
        self.0
    }
}

impl From<u64> for Cost {
    fn from(value: u64) -> Self {
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
