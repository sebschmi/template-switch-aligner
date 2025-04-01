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
    length: usize,
    optional: bool,

    previous: Vec<CharacterKind>,
    current: Option<CharacterKind>,
    next_source: usize,
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

#[derive(Debug, Clone)]
pub enum CharacterKind {
    Source {
        column: SourceColumn,
        kind: SourceCharacterKind,
    },
    Gap {
        identifier: GapIdentifier,
    },
}

#[derive(Debug, Clone)]
pub enum SourceCharacterKind {
    Source { optional: bool },
    Copy { depth: usize },
    Skipped,
}

pub enum AlignedCharacter {
    Source { coordinates: SourceCoordinates },
    Gap { identifier: GapIdentifier },
    Blank,
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
                    Character::new_source(column, data, false)
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
                            SourceCharacterKind::Source { .. } | SourceCharacterKind::Skipped => {
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
        let length = rows[row].length();
        assert!(length > 0);

        let mut result = Self {
            row,
            length,
            optional: rows[row].optional(),

            previous: Default::default(),
            current: Default::default(),
            next_source: 0,
            active_copies: Default::default(),
        };
        result.advance();
        result
    }

    pub fn push_suffix_copy(&mut self, length: usize) {
        assert!(length > 0);
        self.active_copies.push(SuffixCopy::new(length));
    }

    pub fn current(&self) -> Option<CharacterKind> {
        if self.length == self.previous.len() {
            assert!(self.active_copies.is_empty());
            return None;
        }

        Some(self.current.clone().unwrap())
    }

    pub fn advance(&mut self) {
        if let Some(current) = self.current.take() {
            self.previous.push(current);
        } else {
            assert_eq!(self.next_source, 0);
        }

        if self.previous.len() < self.length {
            if let Some(active_copy) = self.active_copies.last_mut() {
                let index = active_copy.apply(self.previous.len());
                if !active_copy.has_remaining_length() {
                    self.active_copies.pop().unwrap();
                }
                self.current = Some(self.previous[index].copy_character())
            } else {
                self.current = Some(CharacterKind::Source {
                    column: self.next_source.into(),
                    kind: SourceCharacterKind::Source {
                        optional: self.optional,
                    },
                });
                self.next_source += 1;
            }
        } else {
            assert!(self.active_copies.is_empty());
        }
    }
}

impl SuffixCopy {
    pub fn new(length: usize) -> Self {
        Self {
            total_length: length,
            remaining_length: length,
        }
    }

    pub fn apply(&mut self, index: usize) -> usize {
        self.remaining_length = self.remaining_length.checked_sub(1).unwrap();
        index.checked_sub(self.total_length).unwrap()
    }

    pub fn has_remaining_length(&self) -> bool {
        self.remaining_length > 0
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

    pub fn new_source(column: SourceColumn, data: Data, optional: bool) -> Self {
        Self::new(
            CharacterKind::Source {
                column,
                kind: SourceCharacterKind::Source { optional },
            },
            data,
        )
    }
}

impl CharacterKind {
    pub fn copy_character(&self) -> Self {
        match self {
            CharacterKind::Source { column, kind } => CharacterKind::Source {
                column: *column,
                kind: match kind {
                    SourceCharacterKind::Source { .. } => SourceCharacterKind::Copy { depth: 0 },
                    SourceCharacterKind::Copy { depth } => {
                        SourceCharacterKind::Copy { depth: *depth + 1 }
                    }
                    SourceCharacterKind::Skipped => unreachable!(),
                },
            },
            CharacterKind::Gap { .. } => unreachable!(),
        }
    }
}
