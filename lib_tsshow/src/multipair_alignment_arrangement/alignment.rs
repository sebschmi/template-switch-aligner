use std::ops::Range;

use tagged_vec::TaggedVec;

use super::coordinates::{AlignmentIndex, SourceColumn, SourceRow};

pub struct Alignment {
    alignment: Vec<AlignmentType>,
    sequence_a: SourceRow,
    sequence_b: SourceRow,
    sequence_a_range: Range<SourceColumn>,
    sequence_b_range: Range<SourceColumn>,
}

/// An alignment type between two sequences `A` and `B`.
#[derive(Debug, Clone, Copy)]
pub enum AlignmentType {
    /// A pair of matching characters.
    Match,
    /// A pair of characters that do not match.
    Substitution,

    /// A gap in sequence `A`.
    ///
    /// The next character in sequence `A` is not consumed, and instead a gap character is inserted.
    GapA,
    /// A gap in sequence `B`.
    ///
    /// The next character in sequence `B` is not consumed, and instead a gap character is inserted.
    GapB,

    /// A skipped character in sequence `A`.
    ///
    /// The next character in sequence `A` is consumed, but not printed.
    SkipA,
    /// A skipped character in sequence `B`.
    ///
    /// The next character in sequence `B` is consumed, but not printed.
    SkipB,
}

pub struct AlignmentIterator {
    index: AlignmentIndex,
    current: usize,
}

impl Alignment {
    pub fn new(
        alignment: Vec<AlignmentType>,
        sequence_a: SourceRow,
        sequence_b: SourceRow,
        sequence_a_range: Range<SourceColumn>,
        sequence_b_range: Range<SourceColumn>,
    ) -> Self {
        Self {
            alignment,
            sequence_a,
            sequence_b,
            sequence_a_range,
            sequence_b_range,
        }
    }

    pub fn new_from_iter<SourceAlignmentType: Into<AlignmentType>>(
        alignment: impl IntoIterator<Item = SourceAlignmentType>,
        sequence_a: SourceRow,
        sequence_b: SourceRow,
        sequence_a_range: Range<SourceColumn>,
        sequence_b_range: Range<SourceColumn>,
    ) -> Self {
        Self::new(
            alignment.into_iter().map(Into::into).collect(),
            sequence_a,
            sequence_b,
            sequence_a_range,
            sequence_b_range,
        )
    }

    pub fn is_between_same_sequence_pair_as(&self, other: &Self) -> bool {
        (self.sequence_a == other.sequence_a && self.sequence_b == other.sequence_b)
            || (self.sequence_a == other.sequence_b && self.sequence_b == other.sequence_a)
    }
}

impl AlignmentIterator {
    pub fn new(index: AlignmentIndex) -> Self {
        Self { index, current: 0 }
    }

    pub fn get<'alignments>(
        &self,
        alignments: &'alignments TaggedVec<AlignmentIndex, Alignment>,
    ) -> Option<AlignmentType> {
        alignments[self.index].alignment.get(self.current).copied()
    }

    pub fn advance(&mut self) {
        self.current += 1;
    }
}
