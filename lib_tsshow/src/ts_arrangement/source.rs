use std::{cmp::Ordering, iter};

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
        copy_depth: Option<usize>,
    },
    Hidden {
        column: SourceColumn,
        copy_depth: Option<usize>,
    },
    Gap {
        copy_depth: Option<usize>,
    },
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
                    result.reference.insert(
                        current_reference_index,
                        SourceChar::Gap {
                            copy_depth: result.query[current_query_index].copy_depth(),
                        },
                    );
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => {
                    result.query.insert(
                        current_query_index,
                        SourceChar::Gap {
                            copy_depth: result.reference[current_reference_index].copy_depth(),
                        },
                    );
                    current_reference_index += 1;
                    current_query_index += 1;
                }
                AlignmentType::PrimarySubstitution | AlignmentType::PrimaryFlankSubstitution => {
                    result.reference[current_reference_index].to_lower_case();
                    result.query[current_query_index].to_lower_case();
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
        let sp2_secondary = match ts_secondary {
            TemplateSwitchSecondary::Reference => {
                self.reference_arrangement_to_source_column(*current_reference_index)
            }
            TemplateSwitchSecondary::Query => {
                self.query_arrangement_to_source_column(*current_query_index)
            }
        } + first_offset;

        let mut sp3_secondary = sp2_secondary;
        let mut primary_inner_length = 0;
        let mut inner_alignment = Alignment::new();

        let anti_primary_gap = loop {
            match alignment.next() {
                Some(AlignmentType::TemplateSwitchExit { anti_primary_gap }) => {
                    break anti_primary_gap;
                }
                Some(alignment_type @ AlignmentType::SecondaryDeletion) => {
                    sp3_secondary -= 1;
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
                    primary_inner_length += 1;
                    inner_alignment.push(alignment_type);
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
                let result = *c;
                c.hide();
                result
            })
            .collect();

        *current_primary_index += primary_inner_length;

        let anti_primary_inner_length = if anti_primary_gap < 0 {
            let duplicate_rev: Vec<_> = anti_primary
                .iter_values()
                .take(current_anti_primary_index.primitive())
                .rev()
                .filter_map(|c| {
                    if c.is_char() {
                        Some(c.make_copy())
                    } else {
                        None
                    }
                })
                .take(usize::try_from(-anti_primary_gap).unwrap())
                .collect();

            anti_primary.splice(
                *current_anti_primary_index..*current_anti_primary_index,
                duplicate_rev.into_iter().rev(),
            );
            0
        } else {
            *current_anti_primary_index += usize::try_from(anti_primary_gap).unwrap();
            usize::try_from(anti_primary_gap).unwrap()
        };

        match primary_inner_length.cmp(&anti_primary_inner_length) {
            Ordering::Less => {
                let delta = anti_primary_inner_length
                    .checked_sub(primary_inner_length)
                    .unwrap();
                primary.splice(
                    *current_primary_index..*current_primary_index,
                    iter::repeat_n(SourceChar::Blank, delta),
                );
                *current_primary_index += delta;
            }
            Ordering::Equal => { /* Do nothing */ }
            Ordering::Greater => {
                let delta = primary_inner_length
                    .checked_sub(anti_primary_inner_length)
                    .unwrap();
                anti_primary.splice(
                    *current_anti_primary_index..*current_anti_primary_index,
                    iter::repeat_n(SourceChar::Blank, delta),
                );
                *current_anti_primary_index += delta;
            }
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
        copy_depth: Option<usize>,
    ) {
        match secondary {
            TemplateSwitchSecondary::Reference => self.insert_reference_gap(column, copy_depth),
            TemplateSwitchSecondary::Query => self.insert_query_gap(column, copy_depth),
        }
    }

    pub fn insert_reference_gap(&mut self, column: ArrangementColumn, copy_depth: Option<usize>) {
        self.reference
            .insert(column, SourceChar::Gap { copy_depth });
        self.query.insert(column, SourceChar::Blank);
    }

    pub fn insert_query_gap(&mut self, column: ArrangementColumn, copy_depth: Option<usize>) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Gap { copy_depth });
    }

    pub fn insert_blank(&mut self, column: ArrangementColumn) {
        self.reference.insert(column, SourceChar::Blank);
        self.query.insert(column, SourceChar::Blank);
    }

    pub fn remove_columns(&mut self, columns: impl IntoIterator<Item = ArrangementColumn> + Clone) {
        self.reference.remove_multi(columns.clone());
        self.query.remove_multi(columns);
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
                SourceChar::Source { column, .. } | SourceChar::Hidden { column, .. }
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

    pub fn reference_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        Self::arrangement_to_source_column(&self.reference, arrangement_column)
    }

    pub fn query_arrangement_to_source_column(
        &self,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        Self::arrangement_to_source_column(&self.query, arrangement_column)
    }

    fn arrangement_to_source_column(
        sequence: &TaggedVec<ArrangementColumn, SourceChar>,
        arrangement_column: ArrangementColumn,
    ) -> SourceColumn {
        // This may also be called on a non-source char.
        assert!(sequence[arrangement_column].is_char());
        sequence
            .iter_values()
            .take(arrangement_column.into())
            .filter(|c| c.is_source_char())
            .count()
            .into()
    }
}

impl SourceChar {
    pub fn new_source(column: impl Into<SourceColumn>) -> Self {
        Self::Source {
            column: column.into(),
            lower_case: false,
            copy_depth: None,
        }
    }

    pub fn is_copy(&self) -> bool {
        match self {
            SourceChar::Source { copy_depth, .. }
            | SourceChar::Hidden { copy_depth, .. }
            | SourceChar::Gap { copy_depth } => copy_depth.is_some(),
            SourceChar::Blank => panic!("Blank has no copy property"),
        }
    }

    pub fn copy_depth(&self) -> Option<usize> {
        match self {
            SourceChar::Source { copy_depth, .. } | SourceChar::Hidden { copy_depth, .. } => {
                *copy_depth
            }
            SourceChar::Gap { copy_depth } => *copy_depth,
            SourceChar::Blank => panic!("Blank has no copy property"),
        }
    }

    pub fn to_lower_case(&mut self) {
        match self {
            Self::Source { lower_case, .. } => *lower_case = true,
            Self::Hidden { .. } | Self::Gap { .. } | Self::Blank => {
                panic!("Not lowercasable: {self:?}")
            }
        }
    }

    pub fn to_upper_case(&mut self) {
        match self {
            Self::Source { lower_case, .. } => *lower_case = false,
            Self::Hidden { .. } | Self::Gap { .. } | Self::Blank => {
                panic!("Not uppercasable: {self:?}")
            }
        }
    }

    pub fn hide(&mut self) {
        match self {
            Self::Source {
                column,
                lower_case,
                copy_depth,
            } => {
                assert!(!*lower_case);
                *self = Self::Hidden {
                    column: *column,
                    copy_depth: *copy_depth,
                };
            }
            Self::Hidden { .. } => unreachable!("Already hidden"),
            Self::Gap { .. } | Self::Blank => unreachable!("Cannot be hidden: {self:?}"),
        }
    }

    pub fn make_copy(&self) -> Self {
        match self {
            Self::Source {
                column,
                lower_case,
                copy_depth,
            } => Self::Source {
                column: *column,
                lower_case: *lower_case,
                copy_depth: Some(copy_depth.map(|copy_depth| copy_depth + 1).unwrap_or(0)),
            },
            Self::Hidden { column, copy_depth } => Self::Hidden {
                column: *column,
                copy_depth: Some(copy_depth.map(|copy_depth| copy_depth + 1).unwrap_or(0)),
            },
            Self::Gap { .. } | Self::Blank => {
                panic!("Should never be copied: {self:?}")
            }
        }
    }
}

impl Char for SourceChar {
    fn source_column(&self) -> SourceColumn {
        match self {
            Self::Source { column, .. } | Self::Hidden { column, .. } => *column,
            Self::Gap { .. } | Self::Blank => panic!("Not a char"),
        }
    }

    fn is_char(&self) -> bool {
        matches!(self, Self::Source { .. } | Self::Hidden { .. })
    }

    fn is_gap(&self) -> bool {
        matches!(self, Self::Gap { .. })
    }

    fn is_blank(&self) -> bool {
        matches!(self, Self::Blank)
    }

    fn is_source_char(&self) -> bool {
        matches!(
            self,
            Self::Source {
                copy_depth: None,
                ..
            } | Self::Hidden {
                copy_depth: None,
                ..
            }
        )
    }

    fn is_hidden(&self) -> bool {
        matches!(self, Self::Hidden { .. })
    }
}
