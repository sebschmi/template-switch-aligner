use std::ops::Add;

/// The cost of an alignment.
///
/// This cost type is not allowed to be negative.
/// This is important for example when using Dijkstra or A* to compute an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
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

impl std::fmt::Display for Cost {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
