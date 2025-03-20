use std::{
    collections::BTreeMap,
    fmt::{Debug, Display},
    iter, mem,
    ops::{Index, IndexMut},
};

use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;
use log::{debug, trace};

#[cfg(test)]
mod tests;

#[derive(Debug)]
pub struct MultipairAlignmentRenderer<SequenceName> {
    sequences: BTreeMap<SequenceName, MultipairAlignmentSequence>,
}

#[derive(Debug)]
pub struct MultipairAlignmentSequence {
    sequence: Vec<Character>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Character {
    Char(char),
    Gap,
    Blank,
}

impl<SequenceName> MultipairAlignmentRenderer<SequenceName> {}

impl<SequenceName: Eq + Ord> MultipairAlignmentRenderer<SequenceName> {
    pub fn new(root_sequence_name: SequenceName, root_sequence: &str) -> Self {
        debug!("Adding root sequence");
        debug!(
            "root_sequence (len: {}): {root_sequence}",
            root_sequence.chars().count()
        );

        Self {
            sequences: [(root_sequence_name, root_sequence.into())]
                .into_iter()
                .collect(),
        }
    }

    pub fn sequence(&self, sequence_name: &SequenceName) -> &MultipairAlignmentSequence {
        self.sequences.get(sequence_name).unwrap()
    }

    /// Append a sequence to the end of an existing rendered sequence.
    ///
    /// Existing gaps and blanks at the end of the existing rendered sequence will not be removed.
    /// The extension will be added after those.
    pub fn extend_sequence(&mut self, sequence_name: &SequenceName, extension: &str) {
        debug!("Extending sequence");
        debug!(
            "extension (len: {}): {extension}",
            extension.chars().count()
        );

        let sequence = self.sequences.get_mut(sequence_name).unwrap();

        // Extend
        sequence.extend_with_string(extension.chars());
        let new_length = sequence.len();

        // Add blanks to other sequences to make all sequences of the same length again
        for (other_sequence_name, other_sequence) in &mut self.sequences {
            if other_sequence_name != sequence_name {
                other_sequence.extend_with_blanks(new_length);
            }
        }
    }

    /// Append a sequence to an existing sequence while aligning it to another sequence.
    #[expect(clippy::too_many_arguments)]
    pub fn extend_sequence_with_alignment(
        &mut self,
        reference_sequence_name: &SequenceName,
        query_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        extension: impl IntoIterator<Item = char>,
        alignment: impl IntoIterator<Item = (usize, AlignmentType)>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        debug!("Extending sequence with alignment");
        debug!("reference_sequence_offset: {reference_sequence_offset}");

        let rendered_sequence_offset = self
            .sequence(reference_sequence_name)
            .translate_alignment_offset(reference_sequence_offset)
            .unwrap_or_else(|| {
                panic!("sequence_offset {reference_sequence_offset} is out of bounds")
            });
        self.sequences
            .get_mut(query_sequence_name)
            .unwrap()
            .prune_blanks(rendered_sequence_offset);
        self.extend_sequence_with_alignment_internal(
            reference_sequence_name,
            query_sequence_name,
            rendered_sequence_offset,
            extension,
            alignment,
            do_lowercasing,
            invert_alignment,
        );
    }

    #[expect(clippy::too_many_arguments)]
    fn extend_sequence_with_alignment_internal(
        &mut self,
        reference_sequence_name: &SequenceName,
        query_sequence_name: &SequenceName,
        rendered_sequence_offset: usize,
        extension: impl IntoIterator<Item = char>,
        alignment: impl IntoIterator<Item = (usize, AlignmentType)>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        let mut reference_gaps = Vec::new();

        let (reference_sequence_name, mut translated_reference_sequence) = self
            .sequences
            .remove_entry(reference_sequence_name)
            .unwrap();
        trace!(
            "translated_reference_sequence.len(): {}",
            translated_reference_sequence.len()
        );
        trace!(
            "translated_reference_sequence_offset: {}",
            rendered_sequence_offset
        );

        let (query_sequence_name, mut translated_query_sequence) =
            self.sequences.remove_entry(query_sequence_name).unwrap();
        assert_eq!(translated_query_sequence.len(), rendered_sequence_offset);

        let mut index = rendered_sequence_offset;
        let mut extension = extension.into_iter();

        for alignment_type in alignment
            .into_iter()
            .flat_map(|(multiplicity, alignment_type)| {
                iter::repeat_n(
                    if invert_alignment {
                        alignment_type.inverted()
                    } else {
                        alignment_type
                    },
                    multiplicity,
                )
            })
        {
            trace!("alignment_type: {alignment_type}");

            while matches!(
                translated_reference_sequence.get(index),
                Some(Character::Blank)
            ) {
                trace!("Skipping blank");
                translated_query_sequence.push(Character::Blank);
                index += 1;
            }

            match alignment_type {
                AlignmentType::PrimaryInsertion
                | AlignmentType::PrimaryFlankInsertion
                | AlignmentType::SecondaryInsertion => {
                    if matches!(
                        translated_reference_sequence.get(index),
                        Some(Character::Gap)
                    ) {
                        index += 1;
                    } else {
                        reference_gaps.push(index);
                    }
                    translated_query_sequence.push(Character::Char(extension.next().unwrap()));
                }
                AlignmentType::PrimaryDeletion
                | AlignmentType::PrimaryFlankDeletion
                | AlignmentType::SecondaryDeletion => {
                    while matches!(
                        translated_reference_sequence.get(index),
                        Some(Character::Gap)
                    ) {
                        translated_query_sequence.push(Character::Blank);
                        index += 1;
                    }
                    translated_query_sequence.push(Character::Gap);
                    index += 1;
                }
                AlignmentType::PrimarySubstitution
                | AlignmentType::PrimaryFlankSubstitution
                | AlignmentType::SecondarySubstitution => {
                    while matches!(
                        translated_reference_sequence.get(index),
                        Some(Character::Gap)
                    ) {
                        translated_query_sequence.push(Character::Blank);
                        index += 1;
                    }
                    let mut extension_character = extension.next().unwrap();

                    if do_lowercasing {
                        extension_character = extension_character.to_ascii_lowercase();
                        if let Character::Char(reference_character) =
                            &mut translated_reference_sequence[index]
                        {
                            *reference_character = reference_character.to_ascii_lowercase();
                        }
                    }

                    translated_query_sequence.push(Character::Char(extension_character));
                    index += 1;
                }
                AlignmentType::PrimaryMatch
                | AlignmentType::PrimaryFlankMatch
                | AlignmentType::SecondaryMatch => {
                    while matches!(
                        translated_reference_sequence.get(index),
                        Some(Character::Gap)
                    ) {
                        trace!("Skipping blank before match");
                        translated_query_sequence.push(Character::Blank);
                        index += 1;
                    }
                    translated_query_sequence.push(Character::Char(extension.next().unwrap()));
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

            assert!(index <= translated_reference_sequence.len());
        }

        assert!(extension.next().is_none());

        translated_query_sequence.extend_with_blanks(translated_reference_sequence.len());
        translated_reference_sequence.insert_gaps(reference_gaps.iter().copied());

        for sequence in self.sequences.values_mut() {
            sequence.insert_blanks(reference_gaps.iter().copied());
        }

        self.sequences
            .insert(reference_sequence_name, translated_reference_sequence);
        self.sequences
            .insert(query_sequence_name, translated_query_sequence);
    }
}

impl<SequenceName: Eq + Ord + Clone> MultipairAlignmentRenderer<SequenceName> {
    #[expect(clippy::too_many_arguments)]
    pub fn add_aligned_sequence(
        &mut self,
        reference_sequence_name: &SequenceName,
        reference_sequence_offset: usize,
        query_sequence_name: SequenceName,
        query_sequence: &str,
        alignment: impl IntoIterator<Item = (usize, AlignmentType)>,
        do_lowercasing: bool,
        invert_alignment: bool,
    ) {
        debug!("Adding aligned sequence");
        debug!("reference_offset: {reference_sequence_offset}");
        debug!(
            "query_sequence (len: {}): {query_sequence}",
            query_sequence.chars().count()
        );
        debug!("invert_alignment: {invert_alignment}");

        assert!(!self.sequences.contains_key(&query_sequence_name));

        let reference_sequence = self.sequences.get_mut(reference_sequence_name).unwrap();
        let index = reference_sequence
            .translate_alignment_offset(reference_sequence_offset)
            .unwrap_or_else(|| {
                panic!("reference_sequence_offset {reference_sequence_offset} is out of bounds")
            });
        self.sequences.insert(
            query_sequence_name.clone(),
            vec![Character::Blank; index].into(),
        );

        self.extend_sequence_with_alignment_internal(
            reference_sequence_name,
            &query_sequence_name,
            index,
            query_sequence.chars(),
            alignment,
            do_lowercasing,
            invert_alignment,
        );
    }

    pub fn add_independent_sequence(&mut self, sequence_name: SequenceName, sequence: &str) {
        assert!(!self.sequences.contains_key(&sequence_name));
        self.sequences.insert(sequence_name, sequence.into());
    }
}

impl<SequenceName: Eq + Ord + Display> MultipairAlignmentRenderer<SequenceName> {
    pub fn render<'name>(
        &self,
        mut output: impl std::io::Write,
        names: impl IntoIterator<Item = &'name SequenceName>,
    ) -> Result<(), std::io::Error>
    where
        SequenceName: 'name,
    {
        let names: Vec<_> = names.into_iter().collect();
        let max_name_len = names
            .iter()
            .map(ToString::to_string)
            .map(|name| name.chars().count())
            .max()
            .unwrap();

        for name in names {
            let sequence = self.sequences.get(name).unwrap();

            let name = name.to_string();
            write!(output, "{name}: ")?;
            for _ in name.len()..max_name_len {
                write!(output, " ")?;
            }

            writeln!(output, "{sequence}")?;
        }

        Ok(())
    }
}

impl<SequenceName: Eq + Ord> MultipairAlignmentRenderer<SequenceName> {
    #[allow(unused)]
    pub fn render_without_names<'name>(
        &self,
        mut output: impl std::io::Write,
        names: impl IntoIterator<Item = &'name SequenceName>,
    ) -> Result<(), std::io::Error>
    where
        SequenceName: 'name,
    {
        let names: Vec<_> = names.into_iter().collect();

        for name in names {
            let sequence = self.sequences.get(name).unwrap();
            writeln!(output, "{sequence}")?;
        }

        Ok(())
    }
}

impl MultipairAlignmentSequence {
    /// Returns the smallest index that skips the first `offset` characters.
    ///
    /// Returns **`None`** if there are less than `offset` characters.
    ///
    /// # Example
    ///
    /// ```rust
    /// let sequence = MultipairAlignmentSequence::from(vec![Character::Blank, Character::Character('A'), Character::Gap, Character::Character('C'), Character::Gap]);
    /// assert_eq!(sequence.translate_alignment_offset(0), 0);
    /// assert_eq!(sequence.translate_alignment_offset(1), 2);
    /// assert_eq!(sequence.translate_alignment_offset(2), 4);
    /// ```
    pub fn translate_alignment_offset(&self, offset: usize) -> Option<usize> {
        if offset == 0 {
            Some(0)
        } else {
            self.sequence
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, character)| matches!(character, Character::Char(_)))
                .nth(offset - 1)
                .map(|(index, _)| index + 1)
        }
    }

    /// Returns the largest index that skips only the first `offset` characters (but not the first `offset + 1` characters).
    ///
    /// Returns **`None`** if there are less than `offset` characters.
    ///
    /// # Example
    ///
    /// ```rust
    /// let sequence = MultipairAlignmentSequence::from(vec![Character::Blank, Character::Character('A'), Character::Gap, Character::Character('C'), Character::Gap]);
    /// assert_eq!(sequence.translate_alignment_offset(0), 1);
    /// assert_eq!(sequence.translate_alignment_offset(1), 3);
    /// assert_eq!(sequence.translate_alignment_offset(2), 5);
    /// ```
    #[expect(unused)]
    pub fn translate_extension_offset(&self, offset: usize) -> Option<usize> {
        self.translate_alignment_offset(offset).map(|offset| {
            self.sequence
                .iter()
                .copied()
                .enumerate()
                .skip(offset)
                .take_while(|(_, character)| !matches!(character, Character::Char(_)))
                .map(|(index, _)| index)
                .last()
                .unwrap_or(offset)
        })
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn characters(&self) -> impl Iterator<Item = char> {
        self.sequence.iter().map(Character::as_char)
    }

    /// Removes blanks from the back of the sequence until the desired length is reached.
    ///
    /// If the desired length is greater than or equal to the current length, then nothing happens.
    /// Panics if any removed character is not a blank.
    pub fn prune_blanks(&mut self, desired_length: usize) {
        while self.len() > desired_length {
            assert_eq!(self.sequence.pop(), Some(Character::Blank));
        }
    }

    /// Adds blanks to the back of the sequence until the desired length is reached.
    ///
    /// If the desired length is lower than or equal to the current length, then nothing happens.
    pub fn extend_with_blanks(&mut self, desired_length: usize) {
        while self.len() < desired_length {
            self.sequence.push(Character::Blank);
        }
    }

    /// Adds the given string to the back of the sequence.
    pub fn extend_with_string(&mut self, string: impl IntoIterator<Item = char>) {
        for character in string {
            self.sequence.push(Character::Char(character));
        }
    }

    /// Adds the given character to the back of the sequence.
    pub fn push(&mut self, character: Character) {
        self.sequence.push(character);
    }

    pub fn get(&self, index: usize) -> Option<&Character> {
        self.sequence.get(index)
    }

    pub fn insert_gaps(&mut self, gaps: impl IntoIterator<Item = usize>) {
        self.multi_insert(Character::Gap, gaps);
    }

    pub fn insert_blanks(&mut self, blanks: impl IntoIterator<Item = usize>) {
        self.multi_insert(Character::Blank, blanks);
    }

    pub fn multi_insert(
        &mut self,
        character: Character,
        positions: impl IntoIterator<Item = usize>,
    ) {
        let original_sequence = mem::take(&mut self.sequence);
        let original_sequence_len = original_sequence.len();

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

        for position in positions {
            debug_assert_eq!(position, original_sequence_len);
            self.sequence.push(character);
        }
    }
}

impl Character {
    fn as_char(&self) -> char {
        match self {
            Character::Char(character) => *character,
            Character::Gap => '-',
            Character::Blank => ' ',
        }
    }
}

impl Display for MultipairAlignmentSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for character in &self.sequence {
            write!(f, "{}", character.as_char())?;
        }

        Ok(())
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
