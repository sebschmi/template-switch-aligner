use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary,
};

use super::{ArrangementCharColumn, ArrangementColumn, SourceColumn, TemplateSwitch};

pub struct TsSourceArrangement {
    reference: Vec<SourceChar>,
    query: Vec<SourceChar>,
}

#[derive(Debug, Clone, Copy)]
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
    Blank,
}

impl TsSourceArrangement {
    pub fn new(
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
        template_switches_out: &mut impl Extend<TemplateSwitch>,
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
                } => template_switches_out.extend([result.align_ts(
                    primary,
                    secondary,
                    first_offset,
                    &mut alignment,
                    &mut current_reference_index,
                    &mut current_query_index,
                )]),

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
        ts_primary: TemplateSwitchPrimary,
        ts_secondary: TemplateSwitchSecondary,
        _first_offset: isize,
        mut alignment: impl Iterator<Item = AlignmentType>,
        current_reference_index: &mut usize,
        current_query_index: &mut usize,
    ) -> TemplateSwitch {
        let sp1_reference = self
            .reference_arrangement_to_char_arrangement_column((*current_reference_index).into());
        let sp1_query =
            self.query_arrangement_to_char_arrangement_column((*current_query_index).into());

        let mut primary_inner_length = 0;
        let mut inner_alignment = Vec::new();

        let anti_primary_gap = loop {
            match alignment.next() {
                Some(AlignmentType::TemplateSwitchExit { anti_primary_gap }) => {
                    break anti_primary_gap;
                }
                Some(
                    alignment_type @ (AlignmentType::SecondaryDeletion
                    | AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch),
                ) => {
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
                }
                Some(alignment_type @ AlignmentType::SecondaryInsertion) => {
                    inner_alignment.push(alignment_type)
                }
                Some(AlignmentType::SecondaryRoot) => { /* Do nothing */ }
                _ => unreachable!(),
            }
        };
        let primary_inner_length = primary_inner_length;

        let (primary, anti_primary, current_primary_index, current_anti_primary_index) =
            match ts_primary {
                TemplateSwitchPrimary::Reference => (
                    &mut self.reference,
                    &mut self.query,
                    current_reference_index,
                    current_query_index,
                ),
                TemplateSwitchPrimary::Query => (
                    &mut self.query,
                    &mut self.reference,
                    current_query_index,
                    current_reference_index,
                ),
            };

        primary
            .iter_mut()
            .skip(*current_primary_index)
            .take(primary_inner_length)
            .for_each(SourceChar::hide);
        *current_primary_index += primary_inner_length;

        if anti_primary_gap < 0 {
            let duplicate: Vec<_> = anti_primary
                .iter()
                .take(*current_anti_primary_index)
                .skip(
                    ((*current_anti_primary_index as isize) - anti_primary_gap)
                        .try_into()
                        .unwrap(),
                )
                .copied()
                .collect();
            anti_primary.splice(
                *current_anti_primary_index..*current_anti_primary_index,
                duplicate,
            );
        } else {
            *current_anti_primary_index += usize::try_from(anti_primary_gap).unwrap();
        }

        let (current_reference_index, current_query_index) = match ts_primary {
            TemplateSwitchPrimary::Reference => (current_primary_index, current_anti_primary_index),
            TemplateSwitchPrimary::Query => (current_anti_primary_index, current_primary_index),
        };

        let sp4_reference = self
            .reference_arrangement_to_char_arrangement_column((*current_reference_index).into());
        let sp4_query =
            self.query_arrangement_to_char_arrangement_column((*current_query_index).into());

        TemplateSwitch {
            primary: ts_primary,
            secondary: ts_secondary,
            sp1_reference,
            sp1_query,
            sp4_reference,
            sp4_query,
            sp2_secondary: todo!(),
            sp3_secondary: todo!(),
            inner_alignment,
        }
    }

    pub fn reference(&self) -> &[SourceChar] {
        &self.reference
    }

    pub fn query(&self) -> &[SourceChar] {
        &self.query
    }

    pub fn reference_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> ArrangementColumn {
        Self::source_to_arrangement_column(&self.reference, column)
    }

    pub fn query_source_to_arrangement_column(&self, column: SourceColumn) -> ArrangementColumn {
        Self::source_to_arrangement_column(&self.query, column)
    }

    fn source_to_arrangement_column(
        sequence: &[SourceChar],
        source_column: SourceColumn,
    ) -> ArrangementColumn {
        sequence
            .iter()
            .enumerate()
            .filter_map(|(i, c)| match c {
                SourceChar::Source { column, .. } | SourceChar::Hidden { column }
                    if *column == source_column =>
                {
                    Some(ArrangementColumn::from(i))
                }
                _ => None,
            })
            .next()
            .unwrap()
    }

    pub fn reference_arrangement_to_char_arrangement_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        Self::arrangement_to_char_arrangement_column(&self.reference, arrangement_column)
    }

    pub fn query_arrangement_to_char_arrangement_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        Self::arrangement_to_char_arrangement_column(&self.query, arrangement_column)
    }

    fn arrangement_to_char_arrangement_column(
        sequence: &[SourceChar],
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
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
            Self::Hidden { .. } | Self::Gap | Self::Blank => panic!("Not lowercaseable: {self:?}"),
        }
    }

    pub fn hide(&mut self) {
        match self {
            Self::Source { column, lower_case }
            | Self::Copy {
                column, lower_case, ..
            } => {
                assert!(!*lower_case);
                *self = Self::Hidden { column: *column };
            }
            Self::Hidden { .. } => unreachable!("Already hidden"),
            Self::Gap | Self::Blank => unreachable!("Cannot be hidden: {self:?}"),
        }
    }
}
