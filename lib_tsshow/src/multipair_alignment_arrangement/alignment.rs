use std::ops::Range;

use tagged_vec::TaggedVec;

use super::{
    coordinates::{AlignmentIndex, SourceColumn, SourceRow},
    sequence::{CharacterKind, CopiedCharactersIterator},
};

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

    /// The start of a copied suffix in sequence `A`.
    CopyA { length: usize },
    /// The start of a copied suffix in sequence `B`.
    CopyB { length: usize },
}

pub struct AlignmentStep {
    alignment_type: AlignmentType,
    sequence_a: SourceRow,
    sequence_b: SourceRow,
    character_a: CharacterKind,
    character_b: CharacterKind,
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

impl AlignmentStep {
    pub fn new(
        alignment_type: AlignmentType,
        sequence_a: SourceRow,
        sequence_b: SourceRow,
        character_a: CharacterKind,
        character_b: CharacterKind,
    ) -> Self {
        Self {
            alignment_type,
            sequence_a,
            sequence_b,
            character_a,
            character_b,
        }
    }

    pub fn is_dominated_by(&self, other: &Self) -> bool {
        assert!(!matches!(
            self.alignment_type,
            AlignmentType::CopyA { .. } | AlignmentType::CopyB { .. }
        ));

        let (self_is_gap, other_is_gap) = match (
            self.sequence_a == other.sequence_a,
            self.sequence_b == other.sequence_b,
            self.sequence_a == other.sequence_b,
            self.sequence_b == other.sequence_a,
        ) {
            (true, _, true, _) | (_, true, _, true) | (true, _, _, true) | (_, true, true, _) => {
                unreachable!("sequence_a and sequence_b cannot be equal within one alignment step")
            }
            (true, true, false, false) | (false, false, true, true) => unreachable!(
                "there is at most one alignment active for each pair of sequences at any point",
            ),
            (true, false, false, false) => (
                matches!(self.alignment_type, AlignmentType::GapA),
                matches!(other.alignment_type, AlignmentType::GapA),
            ),
            (false, true, false, false) => (
                matches!(self.alignment_type, AlignmentType::GapB),
                matches!(other.alignment_type, AlignmentType::GapB),
            ),
            (false, false, true, false) => (
                matches!(self.alignment_type, AlignmentType::GapA),
                matches!(other.alignment_type, AlignmentType::GapB),
            ),
            (false, false, false, true) => (
                matches!(self.alignment_type, AlignmentType::GapB),
                matches!(other.alignment_type, AlignmentType::GapA),
            ),
            (false, false, false, false) => return false,
        };

        if self_is_gap != other_is_gap {
            other_is_gap
        } else {
            false
        }
    }
}

impl AlignmentIterator {
    pub fn new(index: AlignmentIndex) -> Self {
        Self { index, current: 0 }
    }

    pub fn alignment_index(&self) -> AlignmentIndex {
        self.index
    }

    pub fn get(&self, alignments: &TaggedVec<AlignmentIndex, Alignment>) -> Option<AlignmentType> {
        alignments[self.index].alignment.get(self.current).copied()
    }

    pub fn get_alignment_step(
        &self,
        alignments: &TaggedVec<AlignmentIndex, Alignment>,
        next_characters: &TaggedVec<SourceRow, CopiedCharactersIterator>,
    ) -> Option<AlignmentStep> {
        let alignment = &alignments[self.index];
        self.get(alignments).map(|alignment_type| {
            AlignmentStep::new(
                alignment_type,
                alignment.sequence_a,
                alignment.sequence_b,
                next_characters[alignment.sequence_a].current().unwrap(),
                next_characters[alignment.sequence_b].current().unwrap(),
            )
        })
    }

    pub fn advance(&mut self) {
        self.current += 1;
    }
}
