use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
};
use tagged_vec::TaggedVec;

use super::{
    character::Char,
    index_types::{ArrangementCharColumn, ArrangementColumn, SourceColumn},
    template_switch::TemplateSwitch,
};

pub struct TsSourceArrangement {
    reference: TaggedVec<ArrangementColumn, SourceChar>,
    query: TaggedVec<ArrangementColumn, SourceChar>,
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
            reference: FromIterator::from_iter((0..reference_length).map(SourceChar::new_source)),
            query: FromIterator::from_iter((0..query_length).map(SourceChar::new_source)),
        };

        let mut current_reference_index = ArrangementColumn::ZERO;
        let mut current_query_index = ArrangementColumn::ZERO;

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
        first_offset: isize,
        mut alignment: impl Iterator<Item = AlignmentType>,
        current_reference_index: &mut ArrangementColumn,
        current_query_index: &mut ArrangementColumn,
    ) -> TemplateSwitch {
        let sp1_reference =
            self.reference_arrangement_to_char_arrangement_column(*current_reference_index);
        let sp1_query = self.query_arrangement_to_char_arrangement_column(*current_query_index);
        let sp2_secondary = usize::try_from(
            match ts_secondary {
                TemplateSwitchSecondary::Reference => current_reference_index.primitive(),
                TemplateSwitchSecondary::Query => current_query_index.primitive(),
            } as isize
                + first_offset,
        )
        .unwrap()
        .into();

        let mut sp3_secondary = sp2_secondary;
        let mut primary_inner_length = 0;
        let mut inner_alignment = Alignment::new();

        let anti_primary_gap = loop {
            match alignment.next() {
                Some(AlignmentType::TemplateSwitchExit { anti_primary_gap }) => {
                    break anti_primary_gap;
                }
                Some(alignment_type @ AlignmentType::SecondaryDeletion) => {
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
                }
                Some(
                    alignment_type @ (AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch),
                ) => {
                    sp3_secondary -= 1;
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
                }
                Some(alignment_type @ AlignmentType::SecondaryInsertion) => {
                    sp3_secondary -= 1;
                    inner_alignment.push(alignment_type)
                }
                Some(AlignmentType::SecondaryRoot) => { /* Do nothing */ }
                _ => unreachable!(),
            }
        };
        let sp3_secondary = sp3_secondary;
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

        let inner = primary
            .iter_values_mut()
            .skip(current_primary_index.primitive())
            .take(primary_inner_length)
            .map(|c| {
                c.hide();
                *c
            })
            .collect();

        *current_primary_index += primary_inner_length;

        if anti_primary_gap < 0 {
            let duplicate: Vec<_> = anti_primary
                .iter_values()
                .take(current_anti_primary_index.primitive())
                .skip(
                    ((current_anti_primary_index.primitive() as isize) - anti_primary_gap)
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

        let sp4_reference =
            self.reference_arrangement_to_char_arrangement_column(*current_reference_index);
        let sp4_query = self.query_arrangement_to_char_arrangement_column(*current_query_index);

        TemplateSwitch {
            primary: ts_primary,
            secondary: ts_secondary,
            sp1_reference,
            sp1_query,
            sp4_reference,
            sp4_query,
            sp2_secondary,
            sp3_secondary,
            inner,
            inner_alignment,
        }
    }

    pub fn secondary(
        &self,
        secondary: TemplateSwitchSecondary,
    ) -> &TaggedVec<ArrangementColumn, SourceChar> {
        match secondary {
            TemplateSwitchSecondary::Reference => self.reference(),
            TemplateSwitchSecondary::Query => self.query(),
        }
    }

    pub fn reference(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        &self.reference
    }

    pub fn query(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        &self.query
    }

    pub fn insert_secondary_gap(
        &mut self,
        secondary: TemplateSwitchSecondary,
        column: ArrangementColumn,
    ) {
        match secondary {
            TemplateSwitchSecondary::Reference => self.insert_reference_gap(column),
            TemplateSwitchSecondary::Query => self.insert_query_gap(column),
        }
    }

    pub fn insert_reference_gap(&mut self, column: ArrangementColumn) {
        self.reference.insert(column, SourceChar::Gap);
        self.query.insert(column, SourceChar::Blank);
    }

    pub fn insert_query_gap(&mut self, column: ArrangementColumn) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Gap);
    }

    pub fn insert_blank(&mut self, column: ArrangementColumn) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Blank);
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
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        source_column: SourceColumn,
    ) -> ArrangementColumn {
        sequence
            .iter_values()
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
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        arrangement_column: ArrangementColumn,
    ) -> ArrangementCharColumn {
        assert!(sequence[arrangement_column].is_char());
        sequence
            .iter_values()
            .take(arrangement_column.into())
            .filter(|c| c.is_char())
            .count()
            .into()
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
            Self::Hidden { .. } | Self::Gap | Self::Blank => panic!("Not lowercasable: {self:?}"),
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

impl Char for SourceChar {
    fn source_column(&self) -> SourceColumn {
        match self {
            Self::Source { column, .. } | Self::Copy { column, .. } | Self::Hidden { column } => {
                *column
            }
            Self::Gap | Self::Blank => panic!("Not a char"),
        }
    }

    fn is_char(&self) -> bool {
        matches!(
            self,
            Self::Source { .. } | Self::Copy { .. } | Self::Hidden { .. }
        )
    }

    fn is_gap(&self) -> bool {
        matches!(self, Self::Gap)
    }

    fn is_blank(&self) -> bool {
        matches!(self, Self::Blank)
    }

    fn is_source_char(&self) -> bool {
        matches!(self, Self::Source { .. } | Self::Hidden { .. })
    }
}
