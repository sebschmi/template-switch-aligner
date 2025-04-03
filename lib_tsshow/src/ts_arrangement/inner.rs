use std::iter;

use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchSecondary,
};
use tagged_vec::TaggedVec;

use crate::ts_arrangement::character::Char;

use super::{
    complement::TsComplementArrangement,
    index_types::{ArrangementColumn, SourceColumn},
    source::{SourceChar, TsSourceArrangement},
    template_switch::TemplateSwitch,
};

pub struct TsInnerArrangement {
    reference_inners: Vec<TaggedVec<ArrangementColumn, InnerChar>>,
    query_inners: Vec<TaggedVec<ArrangementColumn, InnerChar>>,
}

#[derive(Debug, Clone, Copy)]
pub enum InnerChar {
    Inner {
        column: SourceColumn,
        lower_case: bool,
    },
    Gap,
    Blank,
}

impl TsInnerArrangement {
    pub fn new(
        source_arrangement: &mut TsSourceArrangement,
        complement_arrangement: &mut TsComplementArrangement,
        template_switches: &[TemplateSwitch],
    ) -> Self {
        let mut result = Self {
            reference_inners: Default::default(),
            query_inners: Default::default(),
        };

        for TemplateSwitch {
            secondary,
            sp2_secondary,
            sp3_secondary,
            inner: source_inner,
            inner_alignment,
            ..
        } in template_switches
        {
            let (mut sp2_secondary, sp3_secondary) = match secondary {
                TemplateSwitchSecondary::Reference => (
                    source_arrangement.reference_source_to_arrangement_column(*sp2_secondary),
                    source_arrangement.reference_source_to_arrangement_column(*sp3_secondary),
                ),
                TemplateSwitchSecondary::Query => (
                    source_arrangement.query_source_to_arrangement_column(*sp2_secondary),
                    source_arrangement.query_source_to_arrangement_column(*sp3_secondary),
                ),
            };

            let mut source_inner = source_inner.iter().rev().copied();
            let mut inner = TaggedVec::<ArrangementColumn, _>::default();
            inner.extend(iter::repeat_n(InnerChar::Blank, sp3_secondary.into()));

            let mut current_arrangement_column = sp3_secondary;
            for alignment_type in inner_alignment.iter_flat_cloned().rev() {
                match alignment_type {
                    AlignmentType::SecondaryInsertion => {
                        while !complement_arrangement.secondary_complement(*secondary)
                            [current_arrangement_column]
                            .is_source_char()
                        {
                            inner.push(InnerChar::Blank);
                            current_arrangement_column += 1;
                        }

                        complement_arrangement
                            .show_secondary_character(*secondary, current_arrangement_column);
                        inner.push(InnerChar::Gap);
                        current_arrangement_column += 1;
                    }
                    AlignmentType::SecondaryDeletion => {
                        loop {
                            let c = complement_arrangement.secondary_complement(*secondary)
                                [current_arrangement_column];

                            if c.is_gap() || c.is_source_char() {
                                break;
                            }

                            inner.push(InnerChar::Blank);
                            current_arrangement_column += 1;
                        }

                        if !complement_arrangement.secondary_complement(*secondary)
                            [current_arrangement_column]
                            .is_gap()
                        {
                            complement_arrangement.insert_secondary_complement_gap(
                                *secondary,
                                current_arrangement_column,
                            );

                            source_arrangement.insert_blank(current_arrangement_column);
                            for existing_inner in result
                                .reference_inners
                                .iter_mut()
                                .chain(&mut result.query_inners)
                            {
                                existing_inner.insert(current_arrangement_column, InnerChar::Blank);
                            }

                            sp2_secondary += 1;
                        }

                        inner.push(source_inner.next().unwrap().into());
                        current_arrangement_column += 1;
                    }
                    AlignmentType::SecondarySubstitution | AlignmentType::SecondaryMatch => {
                        while !source_arrangement.secondary(*secondary)[current_arrangement_column]
                            .is_source_char()
                        {
                            inner.push(InnerChar::Blank);
                            current_arrangement_column += 1;
                        }

                        complement_arrangement
                            .show_secondary_character(*secondary, current_arrangement_column);

                        let mut inner_char: InnerChar = source_inner.next().unwrap().into();
                        if alignment_type == AlignmentType::SecondarySubstitution {
                            complement_arrangement
                                .secondary_to_lower_case(*secondary, current_arrangement_column);
                            inner_char.to_lower_case();
                        }

                        inner.push(inner_char);
                        current_arrangement_column += 1;
                    }
                    _ => unreachable!(),
                }
            }

            // We skip further secondary non-source chars for the assertion below.
            while !source_arrangement.secondary(*secondary)[current_arrangement_column]
                .is_source_char()
            {
                current_arrangement_column += 1;
            }
            assert_eq!(current_arrangement_column, sp2_secondary);

            let suffix_blanks =
                iter::repeat_n(InnerChar::Blank, source_arrangement.reference().len())
                    .skip(inner.len());
            inner.extend(suffix_blanks);

            match secondary {
                TemplateSwitchSecondary::Reference => result.reference_inners.push(inner),
                TemplateSwitchSecondary::Query => result.query_inners.push(inner),
            }
        }

        result
    }
}

impl InnerChar {
    pub fn to_lower_case(&mut self) {
        match self {
            InnerChar::Inner { lower_case, .. } => *lower_case = true,
            InnerChar::Gap | InnerChar::Blank => panic!("Not lowercasable"),
        }
    }
}

impl From<SourceChar> for InnerChar {
    fn from(value: SourceChar) -> Self {
        match value {
            SourceChar::Source { column, lower_case } => Self::Inner { column, lower_case },
            SourceChar::Copy { .. } | SourceChar::Hidden { .. } => {
                panic!("Cannot be translated into InnerChar")
            }
            SourceChar::Gap => Self::Gap,
            SourceChar::Blank => Self::Blank,
        }
    }
}
