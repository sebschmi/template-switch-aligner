use std::{fmt::Display, ops::Range};

use serde::Deserialize;

#[derive(Debug, Clone, Eq, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct AlignmentRange {
    offset: AlignmentCoordinates,
    limit: AlignmentCoordinates,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub struct AlignmentCoordinates {
    reference: usize,
    query: usize,
}

impl AlignmentRange {
    pub fn new_complete(reference: usize, query: usize) -> Self {
        Self {
            offset: AlignmentCoordinates::new_zero(),
            limit: AlignmentCoordinates::new(reference, query),
        }
    }

    pub fn new_offset_limit(offset: AlignmentCoordinates, limit: AlignmentCoordinates) -> Self {
        Self { offset, limit }
    }

    pub fn reference_offset(&self) -> usize {
        self.offset.reference
    }

    pub fn query_offset(&self) -> usize {
        self.offset.query
    }

    pub fn reference_limit(&self) -> usize {
        self.limit.reference
    }

    pub fn query_limit(&self) -> usize {
        self.limit.query
    }

    pub fn reference_range(&self) -> Range<usize> {
        self.offset.reference..self.limit.reference
    }

    pub fn query_range(&self) -> Range<usize> {
        self.offset.query..self.limit.query
    }

    pub fn move_offsets_left(&self) -> Self {
        Self::new_offset_limit(
            AlignmentCoordinates::new(
                self.offset.reference.checked_sub(1).unwrap(),
                self.offset.query.checked_sub(1).unwrap(),
            ),
            self.limit,
        )
    }

    pub fn move_offsets_right(&self) -> Self {
        assert!(self.offset.reference < self.limit.reference);
        assert!(self.offset.query < self.limit.query);

        Self::new_offset_limit(
            AlignmentCoordinates::new(self.offset.reference + 1, self.offset.query + 1),
            self.limit,
        )
    }

    pub fn move_limits_left(&self) -> Self {
        assert!(self.offset.reference < self.limit.reference);
        assert!(self.offset.query < self.limit.query);

        Self::new_offset_limit(
            self.offset,
            AlignmentCoordinates::new(self.limit.reference - 1, self.limit.query - 1),
        )
    }

    pub fn move_limits_right(&self) -> Self {
        Self::new_offset_limit(
            self.offset,
            AlignmentCoordinates::new(self.limit.reference + 1, self.limit.query + 1),
        )
    }
}

impl AlignmentCoordinates {
    pub fn new(reference: usize, query: usize) -> Self {
        Self { reference, query }
    }

    pub fn new_zero() -> Self {
        Self::new(0, 0)
    }
}

impl Display for AlignmentRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "R: {}..{}; Q: {}..{}",
            self.offset.reference, self.limit.reference, self.offset.query, self.limit.query
        )
    }
}
