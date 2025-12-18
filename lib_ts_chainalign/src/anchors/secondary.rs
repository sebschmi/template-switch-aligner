use std::fmt::Display;

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates,
        ts_kind::{TsDescendant, TsKind},
    },
    anchors::primary::PrimaryAnchor,
};

/// A secondary anchor.
///
/// This is an anchor between the ancestor in reverse direction and the descendant in forward direction.
///
/// The anchor is ordered by its minimum ordinate first, then by its ancestor ordinate and finally by its descendant ordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SecondaryAnchor {
    /// Ancestor right index in the forward sequence.
    pub(super) ancestor: usize,
    /// Descendant left index in the forward sequence.
    pub(super) descendant: usize,
}

impl SecondaryAnchor {
    pub fn new(ancestor: usize, descendant: usize) -> Self {
        Self {
            ancestor,
            descendant,
        }
    }

    pub fn new_from_start(alignment_coordinates: &AlignmentCoordinates) -> Self {
        Self::new(
            alignment_coordinates.secondary_ordinate_ancestor().unwrap(),
            alignment_coordinates
                .secondary_ordinate_descendant()
                .unwrap(),
        )
    }

    pub fn start(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        AlignmentCoordinates::Secondary {
            ancestor: self.ancestor,
            descendant: self.descendant,
            ts_kind,
        }
    }

    pub fn end(&self, ts_kind: TsKind, k: usize) -> AlignmentCoordinates {
        AlignmentCoordinates::Secondary {
            ancestor: self.ancestor.checked_sub(k).unwrap(),
            descendant: self.descendant + k,
            ts_kind,
        }
    }

    pub fn chaining_gaps(
        &self,
        second: &Self,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<(usize, usize)> {
        let gap_start = self.end(ts_kind, k);
        let gap_end = second.start(ts_kind);

        let gap1 = gap_start
            .secondary_ordinate_ancestor()
            .unwrap()
            .checked_sub(gap_end.secondary_ordinate_ancestor().unwrap())?;
        let gap2 = gap_end
            .secondary_ordinate_descendant()
            .unwrap()
            .checked_sub(gap_start.secondary_ordinate_descendant().unwrap())?;

        Some((gap1, gap2))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &PrimaryAnchor,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(ts_kind, k);
        let gap_end = second.start();

        let gap_start = gap_start.secondary_ordinate_descendant().unwrap();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_end.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => gap_end.primary_ordinate_b().unwrap(),
        };

        gap_end.checked_sub(gap_start)
    }

    pub fn chaining_jump_gap_from_start(
        &self,
        start: AlignmentCoordinates,
        ts_kind: TsKind,
    ) -> usize {
        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => start.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => start.primary_ordinate_b().unwrap(),
        };
        let gap_end = self.start(ts_kind).secondary_ordinate_descendant().unwrap();

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn chaining_jump_gap_to_end(
        &self,
        end: AlignmentCoordinates,
        ts_kind: TsKind,
        k: usize,
    ) -> usize {
        let gap_start = self
            .end(ts_kind, k)
            .secondary_ordinate_descendant()
            .unwrap();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => end.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => end.primary_ordinate_b().unwrap(),
        };

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.ancestor - 1 == successor.ancestor && self.descendant + 1 == successor.descendant
    }

    /// Returns the length of the 2-3 alignment of a TS that starts in `self` and ends in `until`.
    ///
    /// The length is the maximum of the difference of the two sequences.
    pub fn ts_length_until(&self, until: &Self, ts_kind: TsKind, k: usize) -> usize {
        let start = self.start(ts_kind);
        let end = until.end(ts_kind, k);

        (start
            .secondary_ordinate_ancestor()
            .unwrap()
            .checked_sub(end.secondary_ordinate_ancestor().unwrap())
            .unwrap())
        .max(
            end.secondary_ordinate_descendant()
                .unwrap()
                .checked_sub(start.secondary_ordinate_descendant().unwrap())
                .unwrap(),
        )
    }
}

impl Display for SecondaryAnchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.ancestor, self.descendant)
    }
}

impl From<(usize, usize)> for SecondaryAnchor {
    fn from(value: (usize, usize)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl Ord for SecondaryAnchor {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ancestor
            .min(self.descendant)
            .cmp(&other.ancestor.min(other.descendant))
            .then_with(|| self.ancestor.cmp(&other.ancestor))
            .then_with(|| self.descendant.cmp(&other.descendant))
    }
}

impl PartialOrd for SecondaryAnchor {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
