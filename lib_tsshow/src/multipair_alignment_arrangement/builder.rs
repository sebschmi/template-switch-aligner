use std::{collections::HashMap, ops::Range};

use tagged_vec::TaggedVec;

use super::{
    alignment::{Alignment, AlignmentIterator, AlignmentType},
    coordinates::{AlignmentIndex, SourceColumn, SourceRow},
    sequence::CopiedCharactersIterator,
};

#[derive(Default)]
pub struct MultipairAlignmentImplementation {
    rows: TaggedVec<SourceRow, Row>,
    alignments: TaggedVec<AlignmentIndex, Alignment>,
}

pub struct Row {
    length: usize,
    alignment_offsets: HashMap<SourceColumn, Vec<AlignmentIndex>>,
}

pub struct ArrangedColumnIterator<'a> {
    multipair_alignment: &'a MultipairAlignmentImplementation,
    active_alignments: Vec<AlignmentIterator>,
    next_characters: TaggedVec<SourceRow, CopiedCharactersIterator>,
}

impl MultipairAlignmentImplementation {
    pub fn new() -> Self {
        Default::default()
    }

    pub fn add_source_row(&mut self, length: usize) -> SourceRow {
        self.rows.push(Row::new(length))
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
    pub fn new(length: usize) -> Self {
        Self {
            length,
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
}

impl<'a> ArrangedColumnIterator<'a> {
    pub fn new(multipair_alignment: &'a MultipairAlignmentImplementation) -> Self {
        Self {
            multipair_alignment,
            active_alignments: Default::default(),
            next_characters: multipair_alignment
                .rows
                .iter_indices()
                .map(|row| CopiedCharactersIterator::new(row, &multipair_alignment.rows))
                .collect(),
        }
    }
}
