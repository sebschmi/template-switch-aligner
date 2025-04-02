use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;

use super::SourceColumn;

pub struct TsSourceArrangement {
    reference: Vec<SourceChar>,
    query: Vec<SourceChar>,
}

#[derive(Debug)]
pub enum SourceChar {
    Source {
        column: SourceColumn,
        lower_case: bool,
    },
    Copy {
        column: SourceColumn,
        depth: usize,
        lower_case: bool,
    },
    Hidden {
        column: SourceColumn,
    },
    Gap,
}

impl TsSourceArrangement {
    pub fn new(
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
    ) -> Self {
        let mut alignment = alignment.into_iter();
        let mut result = Self {
            reference: Vec::from_iter((0..reference_length).map(SourceChar::new_source)),
            query: Vec::from_iter((0..query_length).map(SourceChar::new_source)),
        };

        let mut current_reference_index = 0;
        let mut current_query_index = 0;

        while let Some(alignment_type) = alignment.next() {
            match alignment_type {
                AlignmentType::PrimaryInsertion | AlignmentType::PrimaryFlankInsertion => {
                    result
                        .reference
                        .insert(current_reference_index, SourceChar::Gap);
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => {
                    result.query.insert(current_query_index, SourceChar::Gap);
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimarySubstitution | AlignmentType::PrimaryFlankSubstitution => {
                    result.reference[current_reference_index].set_lower_case();
                    result.query[current_query_index].set_lower_case();
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimaryMatch | AlignmentType::PrimaryFlankMatch => {
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::TemplateSwitchEntrance {
                    primary,
                    secondary,
                    first_offset,
                } => result.align_ts(
                    &mut alignment,
                    &mut current_reference_index,
                    &mut current_query_index,
                ),

                AlignmentType::Root | AlignmentType::PrimaryReentry => { /* Do nothing */ }
                AlignmentType::TemplateSwitchExit { .. }
                | AlignmentType::SecondaryInsertion
                | AlignmentType::SecondaryDeletion
                | AlignmentType::SecondarySubstitution
                | AlignmentType::SecondaryMatch
                | AlignmentType::SecondaryRoot
                | AlignmentType::PrimaryShortcut { .. } => unreachable!(),
            }
        }

        result
    }

    fn align_ts(
        &mut self,
        mut alignment: impl Iterator<Item = AlignmentType>,
        current_reference_index: &mut usize,
        current_query_index: &mut usize,
    ) {
        todo!()
    }
}

impl SourceChar {
    pub fn new_source(column: impl Into<SourceColumn>) -> Self {
        Self::Source {
            column: column.into(),
            lower_case: false,
        }
    }

    pub fn set_lower_case(&mut self) {
        match self {
            Self::Source { lower_case, .. } | Self::Copy { lower_case, .. } => *lower_case = true,
            Self::Hidden { .. } | Self::Gap => panic!("Not lowercaseable: {self:?}"),
        }
    }
}
