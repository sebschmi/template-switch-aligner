use std::ops::Add;

/// The score of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Score(i64);

impl Score {
    pub const MIN: Score = Score(i64::MIN);
    pub const ZERO: Score = Score(0);

    pub fn as_i64(&self) -> i64 {
        self.0
    }
}

impl From<i64> for Score {
    fn from(value: i64) -> Self {
        Self(value)
    }
}

impl Add for Score {
    type Output = Score;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}
