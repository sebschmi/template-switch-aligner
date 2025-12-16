use std::{
    fmt::Display,
    ops::{Add, Sub},
};

use num_traits::{CheckedSub, Zero};

type IndexType = u32;

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct AnchorIndex(IndexType);

impl AnchorIndex {
    pub fn as_usize(self) -> usize {
        self.0.try_into().unwrap()
    }
}

impl From<IndexType> for AnchorIndex {
    fn from(value: IndexType) -> Self {
        Self(value)
    }
}

impl From<usize> for AnchorIndex {
    fn from(value: usize) -> Self {
        Self(value.try_into().unwrap())
    }
}

impl Add for AnchorIndex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0.checked_add(rhs.0).unwrap())
    }
}

impl Sub for AnchorIndex {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0.checked_sub(rhs.0).unwrap())
    }
}

impl Add<IndexType> for AnchorIndex {
    type Output = Self;

    fn add(self, rhs: IndexType) -> Self::Output {
        self + Self::from(rhs)
    }
}

impl Sub<IndexType> for AnchorIndex {
    type Output = Self;

    fn sub(self, rhs: IndexType) -> Self::Output {
        self - Self::from(rhs)
    }
}

impl CheckedSub for AnchorIndex {
    fn checked_sub(&self, v: &Self) -> Option<Self> {
        self.0.checked_sub(v.0).map(Self)
    }
}

impl Zero for AnchorIndex {
    fn zero() -> Self {
        Self(0)
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl Display for AnchorIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
