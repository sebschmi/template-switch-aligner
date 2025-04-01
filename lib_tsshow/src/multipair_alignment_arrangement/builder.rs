use std::{collections::HashMap, ops::Range};

use tagged_vec::TaggedVec;

use super::{
    alignment::{Alignment, AlignmentIterator, AlignmentType},
    coordinates::{AlignmentIndex, SourceColumn, SourceRow},
    sequence::{CharacterKind, CopiedCharactersIterator},
};

#[derive(Default)]
pub struct MultipairAlignmentImplementation {
    rows: TaggedVec<SourceRow, Row>,
    alignments: TaggedVec<AlignmentIndex, Alignment>,
}

pub struct Row {
    length: usize,
    optional: bool,
    alignment_offsets: HashMap<SourceColumn, Vec<AlignmentIndex>>,
}

pub struct ArrangedColumnIterator<'a> {
    multipair_alignment: &'a MultipairAlignmentImplementation,
    alignment_activated_flags: TaggedVec<AlignmentIndex, bool>,
    active_alignments: Vec<AlignmentIterator>,
    next_characters: TaggedVec<SourceRow, CopiedCharactersIterator>,
}

impl MultipairAlignmentImplementation {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn add_source_row(&mut self, length: usize, optional: bool) -> SourceRow {
        self.rows.push(Row::new(length, optional))
    }

    pub fn add_alignment<SourceAlignmentType: Into<AlignmentType>>(
        &mut self,
        alignment: impl IntoIterator<Item = SourceAlignmentType>,
        sequence_a: SourceRow,
        sequence_b: SourceRow,
        sequence_a_range: Range<SourceColumn>,
        sequence_b_range: Range<SourceColumn>,
    ) {
        debug_assert_ne!(sequence_a, sequence_b);

        let sequence_a_offset = sequence_a_range.start;
        let sequence_b_offset = sequence_b_range.start;
        let alignment = Alignment::new_from_iter(
            alignment,
            sequence_a,
            sequence_b,
            sequence_a_range,
            sequence_b_range,
        );

        for existing_alignment in &self.alignments {
            debug_assert!(!alignment.is_between_same_sequence_pair_as(existing_alignment));
        }

        let alignment = self.alignments.push(alignment);

        self.rows[sequence_a].insert_alignment_offset(sequence_a_offset, alignment);
        self.rows[sequence_a].insert_alignment_offset(sequence_b_offset, alignment);
    }

    pub fn iter_arranged_columns(&self) -> ArrangedColumnIterator {
        ArrangedColumnIterator::new(self)
    }
}

impl Row {
    pub fn new(length: usize, optional: bool) -> Self {
        Self {
            length,
            optional,
            alignment_offsets: Default::default(),
        }
    }

    pub fn insert_alignment_offset(&mut self, column: SourceColumn, alignment: AlignmentIndex) {
        if let Some(alignments) = self.alignment_offsets.get_mut(&column) {
            alignments.push(alignment);
        } else {
            self.alignment_offsets.insert(column, vec![alignment]);
        }
    }
}

impl Row {
    pub fn length(&self) -> usize {
        self.length
    }

    pub fn optional(&self) -> bool {
        self.optional
    }
}

impl<'a> ArrangedColumnIterator<'a> {
    pub fn new(multipair_alignment: &'a MultipairAlignmentImplementation) -> Self {
        let mut alignment_activated_flags: TaggedVec<_, _> =
            vec![false; multipair_alignment.alignments.len()].into();
        let mut active_alignments = Vec::new();

        for row in &multipair_alignment.rows {
            for alignment in row
                .alignment_offsets
                .iter()
                .filter(|(offset, _)| **offset == SourceColumn::from(0))
                .flat_map(|(_, alignments)| alignments.iter().copied())
            {
                if !alignment_activated_flags[alignment] {
                    alignment_activated_flags[alignment] = true;
                    active_alignments.push(AlignmentIterator::new(alignment));

                    todo!(
                        "Only activate alignment when both rows reach it. Assert that both rows are at the correct starting character when reaching."
                    );
                }
            }
        }

        Self {
            multipair_alignment,
            alignment_activated_flags,
            active_alignments,
            next_characters: multipair_alignment
                .rows
                .iter_indices()
                .map(|row| CopiedCharactersIterator::new(row, &multipair_alignment.rows))
                .collect(),
        }
    }

    fn find_dominating_alignments(&self) -> Vec<AlignmentIndex> {
        for alignment1 in &self.active_alignments {}

        todo!()
    }
}

impl Iterator for ArrangedColumnIterator<'_> {
    type Item = TaggedVec<SourceRow, CharacterKind>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut advance_characters: TaggedVec<SourceRow, _> =
            vec![true; self.multipair_alignment.rows.len()].into();

        todo!()
    }
}
