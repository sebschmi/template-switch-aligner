/// A heuristic range within which the start and end of the TS can be shifted without increasing cost.
///
/// The range may not be maximal, but is required to be correct.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EqualCostRange {
    /// How much the start of the TS can be shifted to the left without increasing cost.
    ///
    /// This number must not be positive.
    pub min_start: i8,

    /// How much the start of the TS can be shifted to the right without increasing cost.
    ///
    /// This number must not be negative.
    pub max_start: i8,

    /// How much the end of the TS can be shifted to the left without increasing cost.
    ///
    /// This number must not be positive.
    pub min_end: i8,

    /// How much the end of the TS can be shifted to the right without increasing cost.
    ///
    /// This number must not be negative.
    pub max_end: i8,
}

impl EqualCostRange {
    pub fn new_invalid() -> Self {
        Self {
            min_start: 1,
            max_start: -1,
            min_end: 1,
            max_end: -1,
        }
    }

    pub fn is_valid(&self) -> bool {
        self.min_start <= 0 && self.max_start >= 0 && self.min_end <= 0 && self.max_end >= 0
    }
}
