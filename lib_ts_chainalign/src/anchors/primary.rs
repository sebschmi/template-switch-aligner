use std::fmt::Display;

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates,
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::secondary::SecondaryAnchor,
};

/// A primary anchor.
///
/// This is an anchor between the two sequences in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its first ordinate and finally by its second ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimaryAnchor {
    pub(super) seq1: usize,
    pub(super) seq2: usize,
}

impl PrimaryAnchor {
    pub fn new(seq1: usize, seq2: usize) -> Self {
        Self { seq1, seq2 }
    }

    pub fn new_from_start(alignment_coordinates: &AlignmentCoordinates) -> Self {
        Self::new(
            alignment_coordinates.primary_ordinate_a().unwrap(),
            alignment_coordinates.primary_ordinate_b().unwrap(),
        )
    }

    pub fn new_from_end(alignment_coordinates: &AlignmentCoordinates, k: usize) -> Self {
        Self::new(
            alignment_coordinates
                .primary_ordinate_a()
                .unwrap()
                .checked_sub(k)
                .unwrap(),
            alignment_coordinates
                .primary_ordinate_b()
                .unwrap()
                .checked_sub(k)
                .unwrap(),
        )
    }

    pub fn start(&self) -> AlignmentCoordinates {
        AlignmentCoordinates::Primary {
            a: self.seq1,
            b: self.seq2,
        }
    }

    pub fn end(&self, k: usize) -> AlignmentCoordinates {
        AlignmentCoordinates::Primary {
            a: self.seq1 + k,
            b: self.seq2 + k,
        }
    }

    pub fn chaining_gaps(&self, second: &Self, k: usize) -> Option<(usize, usize)> {
        let gap_start = self.end(k);
        let gap_end = second.start();
        primary_chaining_gaps(gap_start, gap_end)
    }

    pub fn chaining_gaps_from_start(&self, start: AlignmentCoordinates) -> (usize, usize) {
        let gap_end = self.start();
        primary_chaining_gaps(start, gap_end)
            .unwrap_or_else(|| panic!("self: {self}, start: {start}"))
    }

    pub fn chaining_gaps_to_end(&self, end: AlignmentCoordinates, k: usize) -> (usize, usize) {
        let gap_start = self.end(k);
        primary_chaining_gaps(gap_start, end)
            .unwrap_or_else(|| panic!("self: {self}, end: {end}, k: {k}"))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &SecondaryAnchor,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(k);
        let gap_end = second.start(ts_kind);

        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_start.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => gap_start.primary_ordinate_b().unwrap(),
        };
        let gap_end = gap_end.secondary_ordinate_descendant().unwrap();

        gap_end.checked_sub(gap_start)
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.seq1 + 1 == successor.seq1 && self.seq2 + 1 == successor.seq2
    }
}

fn primary_chaining_gaps(
    gap_start: AlignmentCoordinates,
    gap_end: AlignmentCoordinates,
) -> Option<(usize, usize)> {
    let gap1 = gap_end
        .primary_ordinate_a()
        .unwrap()
        .checked_sub(gap_start.primary_ordinate_a().unwrap())?;
    let gap2 = gap_end
        .primary_ordinate_b()
        .unwrap()
        .checked_sub(gap_start.primary_ordinate_b().unwrap())?;

    Some((gap1, gap2))
}

impl Display for PrimaryAnchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.seq1, self.seq2)
    }
}

impl From<(usize, usize)> for PrimaryAnchor {
    fn from(value: (usize, usize)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl Ord for PrimaryAnchor {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.seq1
            .min(self.seq2)
            .cmp(&other.seq1.min(other.seq2))
            .then_with(|| self.seq1.cmp(&other.seq1))
            .then_with(|| self.seq2.cmp(&other.seq2))
    }
}

impl PartialOrd for PrimaryAnchor {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
