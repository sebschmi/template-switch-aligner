use std::{
    collections::BTreeMap,
    fmt::Debug,
    iter, mem,
    ops::{Index, IndexMut},
};

use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;

#[derive(Debug)]
pub struct MultipairAlignmentRenderer<SequenceName> {
    sequences: BTreeMap<SequenceName, MultipairAlignmentSequence>,
}

#[derive(Debug)]
struct MultipairAlignmentSequence {
    sequence: Vec<Character>,
}

#[derive(Debug, Clone, Copy)]
enum Character {
    Char(char),
    Gap,
    Blank,
}

impl<SequenceName> MultipairAlignmentRenderer<SequenceName> {}

impl<SequenceName: Eq + Ord> MultipairAlignmentRenderer<SequenceName> {
    pub fn new(root_sequence_name: SequenceName, root_sequence: &str) -> Self {
        Self {
            sequences: [(root_sequence_name, root_sequence.into())]
                .into_iter()
                .collect(),
        }
    }
}

impl<SequenceName: Eq + Ord + Clone> MultipairAlignmentRenderer<SequenceName> {
    #[expect(clippy::too_many_arguments)]
    pub fn add_aligned_sequence(
        &mut self,
        reference_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        new_sequence_name: SequenceName,
        new_sequence: &str,
        alignment: &[(usize, AlignmentType)],
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        assert!(!self.sequences.contains_key(&new_sequence_name));

        let mut new_sequence = new_sequence.chars();
        let reference_sequence = self.sequences.get_mut(reference_sequence_name).unwrap();
        let mut index = reference_sequence.nth_character_index(reference_sequence_offset);
        let mut translated_new_sequence = vec![Character::Blank; index];
        let mut reference_gaps = Vec::new();

        for alignment_type in
            alignment
                .iter()
                .copied()
                .flat_map(|(multiplicity, alignment_type)| {
                    iter::repeat(if invert_alignment {
                        alignment_type.inverted()
                    } else {
                        alignment_type
                    })
                    .take(multiplicity)
                })
        {
            while matches!(reference_sequence[index], Character::Gap | Character::Blank) {
                translated_new_sequence.push(Character::Blank);
                index += 1;
            }

            match alignment_type {
                AlignmentType::PrimaryInsertion
                | AlignmentType::PrimaryFlankInsertion
                | AlignmentType::SecondaryInsertion => {
                    reference_gaps.push(index);
                    translated_new_sequence.push(Character::Char(new_sequence.next().unwrap()));
                }
                AlignmentType::PrimaryDeletion
                | AlignmentType::PrimaryFlankDeletion
                | AlignmentType::SecondaryDeletion => {
                    translated_new_sequence.push(Character::Gap);
                    index += 1;
                }
                AlignmentType::PrimarySubstitution
                | AlignmentType::PrimaryFlankSubstitution
                | AlignmentType::SecondarySubstitution => {
                    let mut new_character = new_sequence.next().unwrap();

                    if do_lowercasing {
                        new_character = new_character.to_ascii_lowercase();
                        if let Character::Char(reference_character) = &mut reference_sequence[index]
                        {
                            *reference_character = reference_character.to_ascii_lowercase();
                        }
                    }

                    translated_new_sequence.push(Character::Char(new_character));
                    index += 1;
                }
                AlignmentType::PrimaryMatch
                | AlignmentType::PrimaryFlankMatch
                | AlignmentType::SecondaryMatch => {
                    translated_new_sequence.push(Character::Char(new_sequence.next().unwrap()));
                    index += 1;
                }
                AlignmentType::Root
                | AlignmentType::PrimaryReentry
                | AlignmentType::TemplateSwitchEntrance { .. }
                | AlignmentType::TemplateSwitchExit { .. }
                | AlignmentType::SecondaryRoot
                | AlignmentType::PrimaryShortcut { .. } => {
                    panic!("Not allowed in rendered alignment: {alignment_type:?}")
                }
            }
        }

        assert!(new_sequence.next().is_none());
        assert!(index <= reference_sequence.len());

        translated_new_sequence
            .extend(iter::repeat(Character::Blank).take(reference_sequence.len() - index));
        reference_sequence.insert_gaps(reference_gaps.iter().copied());

        for (name, sequence) in self.sequences.iter_mut() {
            if name != reference_sequence_name {
                sequence.insert_blanks(reference_gaps.iter().copied());
            }
        }

        self.sequences
            .insert(new_sequence_name, translated_new_sequence.into());
    }
}

impl<SequenceName: Eq + Ord + Debug> MultipairAlignmentRenderer<SequenceName> {
    pub fn render(&self, mut output: impl std::io::Write) -> Result<(), std::io::Error> {
        write!(output, "{self:?}")
    }
}

impl MultipairAlignmentSequence {
    fn nth_character_index(&self, n: usize) -> usize {
        self.sequence
            .iter()
            .enumerate()
            .filter(|(_, character)| matches!(character, Character::Char(_)))
            .nth(n)
            .unwrap()
            .0
    }

    fn len(&self) -> usize {
        self.sequence.len()
    }

    fn insert_gaps(&mut self, gaps: impl IntoIterator<Item = usize>) {
        self.multi_insert(Character::Gap, gaps);
    }

    fn insert_blanks(&mut self, blanks: impl IntoIterator<Item = usize>) {
        self.multi_insert(Character::Blank, blanks);
    }

    fn multi_insert(&mut self, character: Character, positions: impl IntoIterator<Item = usize>) {
        let original_sequence = mem::take(&mut self.sequence);

        let mut positions = positions.into_iter().peekable();
        let original_characters = original_sequence.into_iter().enumerate();

        for (index, original_character) in original_characters {
            while let Some(position) = positions.peek().copied() {
                if position <= index {
                    self.sequence.push(character);
                    positions.next().unwrap();
                } else {
                    break;
                }
            }

            self.sequence.push(original_character);
        }
    }
}

impl From<&str> for MultipairAlignmentSequence {
    fn from(value: &str) -> Self {
        Self {
            sequence: value.chars().map(Character::Char).collect(),
        }
    }
}

impl From<Vec<Character>> for MultipairAlignmentSequence {
    fn from(value: Vec<Character>) -> Self {
        Self { sequence: value }
    }
}

impl Index<usize> for MultipairAlignmentSequence {
    type Output = <Vec<Character> as Index<usize>>::Output;

    fn index(&self, index: usize) -> &Self::Output {
        self.sequence.index(index)
    }
}

impl IndexMut<usize> for MultipairAlignmentSequence {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.sequence.index_mut(index)
    }
}
