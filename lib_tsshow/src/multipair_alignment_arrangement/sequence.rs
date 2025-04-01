use tagged_vec::TaggedVec;

use super::{
    builder::Row,
    coordinates::{GapIdentifier, SourceColumn, SourceCoordinates, SourceRow},
};

pub struct AlignedSequence<Data> {
    sequence: Vec<Character<Data>>,
}

pub struct CopiedCharactersIterator {
    row: SourceRow,
    current: usize,
    length: usize,
    active_copies: Vec<SuffixCopy>,
}

struct SuffixCopy {
    total_length: usize,
    remaining_length: usize,
}

pub struct Character<Data> {
    data: Data,
    kind: CharacterKind,
    aligned_characters: Vec<AlignedCharacter>,
}

pub enum CharacterKind {
    Source {
        column: SourceColumn,
        kind: SourceCharacterKind,
    },
    Gap {
        identifier: GapIdentifier,
    },
}

pub enum SourceCharacterKind {
    Source,
    Copy { depth: usize },
    Skipped,
}

pub enum AlignedCharacter {
    Source { coordinates: SourceCoordinates },
    Gap { identifier: GapIdentifier },
}

impl<Data> AlignedSequence<Data> {
    pub fn new_source_sequence(
        length: usize,
        mut data_generator: impl FnMut(SourceColumn) -> Data,
    ) -> Self {
        Self {
            sequence: (0..length)
                .map(|column| {
                    let column = column.into();
                    let data = data_generator(column);
                    Character::new_source(column, data)
                })
                .collect(),
        }
    }

    /// Duplicates a suffix of a sequence.
    ///
    /// The parameter `length` decides the length of the suffix to duplicate.
    /// The suffix only includes source characters, and skips e.g. gaps, not counting them into the length.
    ///
    /// **Example:** Duplicating the suffix of length `2` within the sequence `AB-A` yields the sequence `AB-ABA`.
    ///
    /// The `data_generator` receives the orignal character that is going to be duplicated along with the [`CharacterKind`] of the duplicate.
    pub fn duplicate_source_suffix(
        &mut self,
        length: usize,
        mut data_generator: impl FnMut(&Character<Data>, &CharacterKind) -> Data,
    ) {
        let extension: Vec<_> = self
            .sequence
            .iter()
            .rev()
            .filter_map(|character| match &character.kind {
                CharacterKind::Source { column, kind } => {
                    let duplicate = CharacterKind::Source {
                        column: *column,
                        kind: match kind {
                            SourceCharacterKind::Source | SourceCharacterKind::Skipped => {
                                SourceCharacterKind::Copy { depth: 0 }
                            }
                            SourceCharacterKind::Copy { depth } => {
                                SourceCharacterKind::Copy { depth: depth + 1 }
                            }
                        },
                    };

                    let data = data_generator(character, &duplicate);
                    Some(Character::new(duplicate, data))
                }
                CharacterKind::Gap { .. } => None,
            })
            .take(length)
            .collect();
        self.sequence.extend(extension);
    }
}

impl CopiedCharactersIterator {
    pub fn new(row: SourceRow, rows: &TaggedVec<SourceRow, Row>) -> Self {
        Self {
            row,
            current: 0,
            length: rows[row].length(),
            active_copies: Default::default(),
        }
    }
}

impl<Data> Character<Data> {
    pub fn new(kind: CharacterKind, data: Data) -> Self {
        Self {
            data,
            kind,
            aligned_characters: Default::default(),
        }
    }

    pub fn new_source(column: SourceColumn, data: Data) -> Self {
        Self::new(
            CharacterKind::Source {
                column,
                kind: SourceCharacterKind::Source,
            },
            data,
        )
    }
}
