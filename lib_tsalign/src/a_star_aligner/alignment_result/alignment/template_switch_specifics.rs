use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::GenomeSequence,
};

use crate::a_star_aligner::{
    alignment_result::{IAlignmentType, alignment::stream::AlignmentStream},
    template_switch_distance::{
        AlignmentType, TemplateSwitchDirection, TemplateSwitchPrimary, TemplateSwitchSecondary,
    },
};

use super::Alignment;

impl Alignment<AlignmentType> {
    /// Moves the start of a template switch one character pair to the left if possible, moving the character pair into the template switch.
    ///
    /// Only works if the alignment preceding the template switch is a match or substitution, and not a flank.
    /// In this case it returns true, and otherwise false.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_start_backwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        mut compact_index: usize,
    ) -> bool {
        let AlignmentType::TemplateSwitchEntrance {
            first_offset,
            primary,
            secondary,
            direction,
            ..
        } = self.alignment[compact_index].1
        else {
            panic!()
        };

        if let Some((_, AlignmentType::PrimaryMatch | AlignmentType::PrimarySubstitution)) =
            self.alignment.get(compact_index - 1)
        {
            // Compute TS inner first indices.
            let mut stream = AlignmentStream::new();
            stream.push_all(self.iter_compact_cloned().take(compact_index));
            let ts_inner_primary_index = match primary {
                TemplateSwitchPrimary::Reference => {
                    stream.head_coordinates().reference() + reference_offset
                }
                TemplateSwitchPrimary::Query => stream.head_coordinates().query() + query_offset,
            };
            let ts_inner_secondary_index = usize::try_from(
                isize::try_from(match secondary {
                    TemplateSwitchSecondary::Reference => {
                        stream.head_coordinates().reference() + reference_offset
                    }
                    TemplateSwitchSecondary::Query => {
                        stream.head_coordinates().query() + query_offset
                    }
                })
                .unwrap()
                .checked_add(first_offset)
                .unwrap(),
            )
            .unwrap();

            // Remove one match or substitution from before the TS.
            let multiplicity = &mut self.alignment[compact_index - 1].0;
            assert!(*multiplicity > 0);
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                compact_index -= 1;
                self.alignment.remove(compact_index);
            }

            // Check if the new inner pair is a match or a substitution.
            let primary_char = primary.get(reference, query)[ts_inner_primary_index - 1].clone();
            let secondary_char = match direction {
                TemplateSwitchDirection::Forward => {
                    secondary.get(reference, query)[ts_inner_secondary_index - 1].clone()
                }
                TemplateSwitchDirection::Reverse => {
                    secondary.get(reference, query)[ts_inner_secondary_index].complement()
                }
            };
            let inner_alignment_type = if primary_char == secondary_char {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            };

            // Insert new inner alignment.
            if self.alignment[compact_index + 1].1 == inner_alignment_type {
                self.alignment[compact_index + 1].0 += 1;
            } else {
                self.alignment
                    .insert(compact_index + 1, (1, inner_alignment_type));
            }

            // If reverse TS, then fix first offset.
            // (If forward TS, the changes to points 1 and 2 cancel out.)
            if direction == TemplateSwitchDirection::Reverse {
                let AlignmentType::TemplateSwitchEntrance { first_offset, .. } =
                    &mut self.alignment[compact_index].1
                else {
                    unreachable!();
                };
                *first_offset += 2;
            }

            // Fix anti-primary gap.
            let Some((_, AlignmentType::TemplateSwitchExit { anti_primary_gap })) = self.alignment
                [compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!();
            };
            *anti_primary_gap += 1;

            true
        } else {
            false
        }
    }

    /// Moves the start of a template switch one character pair to the right if possible, moving the character pair out of the template switch.
    ///
    /// Only works if the first alignment inside the template switch is a match or substitution.
    /// In this case it returns true, and otherwise false.
    ///
    /// Assumes that the template switch is not preceded by a flank.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_start_forwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        mut compact_index: usize,
    ) -> bool {
        let AlignmentType::TemplateSwitchEntrance { direction, .. } =
            self.alignment[compact_index].1
        else {
            panic!()
        };

        // Assert that no flanks are involved.
        assert!(
            self.alignment
                .get(compact_index - 1)
                .map(|(_, alignment_type)| !matches!(
                    alignment_type,
                    AlignmentType::PrimaryFlankDeletion
                        | AlignmentType::PrimaryFlankInsertion
                        | AlignmentType::PrimaryFlankSubstitution
                        | AlignmentType::PrimaryFlankMatch
                ))
                .unwrap_or(true)
        );

        if let Some((_, AlignmentType::SecondaryMatch | AlignmentType::SecondarySubstitution)) =
            self.alignment.get(compact_index + 1)
        {
            // Compute TS outer first indices.
            let mut stream = AlignmentStream::new();
            stream.push_all(self.iter_compact_cloned().take(compact_index));
            let ts_outer_reference_index = stream.head_coordinates().reference() + reference_offset;
            let ts_outer_query_index = stream.head_coordinates().query() + query_offset;

            // Remove one match or substitution from inside the TS.
            let multiplicity = &mut self.alignment[compact_index + 1].0;
            assert!(*multiplicity > 0);
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                self.alignment.remove(compact_index + 1);
            }

            // Check if the new outer pair is a match or a substitution.
            let reference_char = reference[ts_outer_reference_index].clone();
            let query_char = query[ts_outer_query_index].clone();
            let outer_alignment_type = if reference_char == query_char {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert new outer alignment.
            if self.alignment[compact_index - 1].1 == outer_alignment_type {
                self.alignment[compact_index - 1].0 += 1;
            } else {
                self.alignment
                    .insert(compact_index - 1, (1, outer_alignment_type));
                compact_index += 1;
            }

            // If reverse TS, then fix first offset.
            // (If forward TS, the changes to points 1 and 2 cancel out.)
            if direction == TemplateSwitchDirection::Reverse {
                let AlignmentType::TemplateSwitchEntrance { first_offset, .. } =
                    &mut self.alignment[compact_index].1
                else {
                    unreachable!();
                };
                *first_offset -= 2;
            }

            // Fix anti-primary gap.
            let Some((_, AlignmentType::TemplateSwitchExit { anti_primary_gap })) = self.alignment
                [compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!();
            };
            *anti_primary_gap -= 1;

            true
        } else {
            false
        }
    }

    /// Moves the end of a template switch one character pair to the right if possible, moving the character pair into the template switch.
    ///
    /// Only works if the alignment preceding the template switch is a match or substitution, and not a flank.
    /// In this case it returns true, and otherwise false.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_end_forwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: usize,
    ) -> bool {
        let AlignmentType::TemplateSwitchEntrance {
            first_offset,
            primary,
            secondary,
            direction,
            ..
        } = self.alignment[compact_index].1
        else {
            panic!()
        };
        let mut exit_index = compact_index
            + self
                .alignment
                .iter()
                .skip(compact_index)
                .enumerate()
                .find(|(_, (_, alignment_type))| {
                    matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. })
                })
                .unwrap()
                .0;
        let inner_secondary_length = self
            .alignment
            .iter()
            .take(exit_index)
            .skip(compact_index + 1)
            .fold(
                0,
                |secondary_length, (multiplicity, alignment_type)| match alignment_type {
                    AlignmentType::SecondaryInsertion => secondary_length,
                    AlignmentType::SecondaryDeletion
                    | AlignmentType::SecondarySubstitution
                    | AlignmentType::SecondaryMatch => secondary_length + *multiplicity,
                    _ => unreachable!(),
                },
            );
        println!("Inner secondary length: {inner_secondary_length}");

        if let Some((_, AlignmentType::PrimaryMatch | AlignmentType::PrimarySubstitution)) =
            self.alignment.get(exit_index + 1)
        {
            // Compute TS inner last indices.
            let mut stream = AlignmentStream::new();
            stream.push_all(self.iter_compact_cloned().take(compact_index));
            stream.clear();
            stream.push_all(
                self.iter_compact_cloned()
                    .take(exit_index + 1)
                    .skip(compact_index),
            );
            let ts_inner_primary_index = match primary {
                TemplateSwitchPrimary::Reference => {
                    stream.head_coordinates().reference() + reference_offset
                }
                TemplateSwitchPrimary::Query => stream.head_coordinates().query() + query_offset,
            };
            let ts_inner_secondary_index = usize::try_from(
                isize::try_from(match secondary {
                    TemplateSwitchSecondary::Reference => {
                        stream.tail_coordinates().reference() + reference_offset
                    }
                    TemplateSwitchSecondary::Query => {
                        stream.tail_coordinates().query() + query_offset
                    }
                })
                .unwrap()
                .checked_add(first_offset)
                .unwrap(),
            )
            .unwrap();
            let ts_inner_secondary_index = match direction {
                TemplateSwitchDirection::Forward => ts_inner_secondary_index
                    .checked_add(inner_secondary_length)
                    .unwrap(),
                TemplateSwitchDirection::Reverse => ts_inner_secondary_index
                    .checked_sub(inner_secondary_length)
                    .unwrap(),
            };

            // Remove one match or substitution from after the TS.
            let multiplicity = &mut self.alignment[exit_index + 1].0;
            assert!(*multiplicity > 0);
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                self.alignment.remove(exit_index + 1);
            }

            // Check if the new inner pair is a match or a substitution.
            println!("Primary index: {ts_inner_primary_index}");
            let primary_char = primary.get(reference, query)[ts_inner_primary_index].clone();
            let secondary_char = match direction {
                TemplateSwitchDirection::Forward => {
                    secondary.get(reference, query)[ts_inner_secondary_index].clone()
                }
                TemplateSwitchDirection::Reverse => {
                    println!("Secondary index: {}", ts_inner_secondary_index - 1);
                    secondary.get(reference, query)[ts_inner_secondary_index - 1].complement()
                }
            };
            let inner_alignment_type = if primary_char == secondary_char {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            };

            // Insert new inner alignment.
            if self.alignment[exit_index - 1].1 == inner_alignment_type {
                self.alignment[exit_index - 1].0 += 1;
            } else {
                self.alignment.insert(exit_index, (1, inner_alignment_type));
                exit_index += 1;
            }

            // Fix anti-primary gap.
            let AlignmentType::TemplateSwitchExit { anti_primary_gap } =
                &mut self.alignment[exit_index].1
            else {
                unreachable!();
            };
            *anti_primary_gap += 1;

            true
        } else {
            false
        }
    }

    /// Moves the end of a template switch one character pair to the left if possible, moving the character pair out of the template switch.
    ///
    /// Only works if the last alignment inside the template switch is a match or substitution.
    /// In this case it returns true, and otherwise false.
    ///
    /// Assumes that the template switch is not succeeded by a flank.
    ///
    /// Compact index identifies the start of the template switch in terms of the compact representation of the alignment.
    /// See [`Self::iter_compact`].
    pub fn move_template_switch_end_backwards<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_offset: usize,
        query_offset: usize,
        compact_index: usize,
    ) -> bool {
        assert!(matches!(
            self.alignment[compact_index].1,
            AlignmentType::TemplateSwitchEntrance { .. }
        ));
        let mut exit_index = compact_index
            + self
                .alignment
                .iter()
                .skip(compact_index)
                .enumerate()
                .find(|(_, (_, alignment_type))| {
                    matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. })
                })
                .unwrap()
                .0;

        // Assert that no flanks are involved.
        assert!(
            self.alignment
                .get(exit_index + 1)
                .map(|(_, alignment_type)| !matches!(
                    alignment_type,
                    AlignmentType::PrimaryFlankDeletion
                        | AlignmentType::PrimaryFlankInsertion
                        | AlignmentType::PrimaryFlankSubstitution
                        | AlignmentType::PrimaryFlankMatch
                ))
                .unwrap_or(true)
        );

        if let Some((_, AlignmentType::SecondaryMatch | AlignmentType::SecondarySubstitution)) =
            self.alignment.get(exit_index - 1)
        {
            // Compute TS outer first indices.
            let mut stream = AlignmentStream::new();
            stream.push_all(self.iter_compact_cloned().take(exit_index + 1));
            let ts_outer_reference_index = stream.head_coordinates().reference() + reference_offset;
            let ts_outer_query_index = stream.head_coordinates().query() + query_offset;

            // Remove one match or substitution from inside the TS.
            let multiplicity = &mut self.alignment[exit_index - 1].0;
            assert!(*multiplicity > 0);
            *multiplicity -= 1;

            // Remove the alignment entry if it has zero multiplicity.
            if *multiplicity == 0 {
                exit_index -= 1;
                self.alignment.remove(exit_index);
            }

            // Check if the new outer pair is a match or a substitution.
            let reference_char = reference[ts_outer_reference_index].clone();
            let query_char = query[ts_outer_query_index].clone();
            let outer_alignment_type = if reference_char == query_char {
                AlignmentType::PrimaryMatch
            } else {
                AlignmentType::PrimarySubstitution
            };

            // Insert new outer alignment.
            if self.alignment[exit_index + 1].1 == outer_alignment_type {
                self.alignment[exit_index + 1].0 += 1;
            } else {
                self.alignment
                    .insert(exit_index + 1, (1, outer_alignment_type));
            }

            // Fix anti-primary gap.
            let Some((_, AlignmentType::TemplateSwitchExit { anti_primary_gap })) = self.alignment
                [compact_index..]
                .iter_mut()
                .find(|(_, alignment_type)| alignment_type.is_template_switch_exit())
            else {
                unreachable!();
            };
            *anti_primary_gap -= 1;

            true
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
        interface::sequence::{GenomeSequence, OwnedGenomeSequence},
    };

    use crate::a_star_aligner::{
        alignment_result::alignment::Alignment,
        template_switch_distance::{
            AlignmentType, EqualCostRange, TemplateSwitchDirection, TemplateSwitchPrimary,
            TemplateSwitchSecondary,
        },
    };

    static START_REFERENCE: &[u8] = b"AGAGAGCTCTAA";
    static START_QUERY: &[u8] = b"AGAGAGCTTTAA";
    static START_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
        &[
            (6, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -6,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 2,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (5, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -4,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (3, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 3,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (4, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: -2,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 4,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (3, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 0,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (1, AlignmentType::SecondarySubstitution),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 5,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
        &[
            (2, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 2,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondarySubstitution),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 6,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
    ];

    static END_REFERENCE: &[u8] = b"AACTCTAGAGAG";
    static END_QUERY: &[u8] = b"AATTCTAGAGAG";
    static END_ALIGNMENTS: &[&[(usize, AlignmentType)]] = &[
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (2, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 2,
                },
            ),
            (6, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (3, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 3,
                },
            ),
            (5, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 4,
                },
            ),
            (4, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (1, AlignmentType::SecondarySubstitution),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 5,
                },
            ),
            (3, AlignmentType::PrimaryMatch),
        ],
        &[
            (1, AlignmentType::PrimaryMatch),
            (
                1,
                AlignmentType::TemplateSwitchEntrance {
                    first_offset: 10,
                    equal_cost_range: EqualCostRange::new_invalid(),
                    primary: TemplateSwitchPrimary::Reference,
                    secondary: TemplateSwitchSecondary::Query,
                    direction: TemplateSwitchDirection::Reverse,
                },
            ),
            (4, AlignmentType::SecondaryMatch),
            (2, AlignmentType::SecondarySubstitution),
            (
                1,
                AlignmentType::TemplateSwitchExit {
                    anti_primary_gap: 6,
                },
            ),
            (2, AlignmentType::PrimaryMatch),
        ],
    ];

    #[test]
    fn move_template_switch_start_backwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
        let mut alignment = Alignment::from(START_ALIGNMENTS[0].to_vec());

        for expected_alignment in &START_ALIGNMENTS[1..] {
            assert!(alignment.move_template_switch_start_backwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        }
    }

    #[test]
    fn move_template_switch_start_forwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(START_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(START_QUERY).unwrap();
        let mut alignment = Alignment::from(START_ALIGNMENTS.last().unwrap().to_vec());

        for expected_alignment in START_ALIGNMENTS.iter().rev().skip(1) {
            assert!(alignment.move_template_switch_start_forwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                2,
                2,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        }
    }

    #[test]
    fn move_template_switch_end_backwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
        let mut alignment = Alignment::from(END_ALIGNMENTS.last().unwrap().to_vec());

        for expected_alignment in END_ALIGNMENTS.iter().rev().skip(1) {
            assert!(alignment.move_template_switch_end_backwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        }
    }

    #[test]
    fn move_template_switch_end_forwards() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(END_REFERENCE).unwrap();
        let query = VectorGenome::from_slice_u8(END_QUERY).unwrap();
        let mut alignment = Alignment::from(END_ALIGNMENTS[0].to_vec());

        for expected_alignment in &END_ALIGNMENTS[1..] {
            assert!(alignment.move_template_switch_end_forwards(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                1,
                1,
                1
            ));
            assert_eq!(alignment, Alignment::from(expected_alignment.to_vec()));
        }
    }
}
