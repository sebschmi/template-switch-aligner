use anyhow::{Result, anyhow};
use compact_genome::io::fasta::FastaRecord;
use log::debug;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::Path,
};
use utf8_chars::BufReadCharsExt;

pub fn parse_pair_fasta_file(
    path: impl AsRef<Path>,
) -> Result<(FastaRecord<String>, FastaRecord<String>)> {
    let mut result = parse_fasta_file(path, true)?;
    let second = result.remove(1);
    let first = result.remove(0);
    Ok((first, second))
}

pub fn parse_single_fasta_file(path: impl AsRef<Path>) -> Result<FastaRecord<String>> {
    Ok(parse_fasta_file(path, false)?.remove(0))
}

fn parse_fasta_file(path: impl AsRef<Path>, pair: bool) -> Result<Vec<FastaRecord<String>>> {
    let path = path.as_ref();
    debug!("Parsing fasta file {path:?}");

    enum State {
        FileStart,
        ParseId,
        ParseComment,
        ParseSequence,
    }

    let mut input = CharacterIterator::new(BufReader::new(
        File::open(path).map_err(|error| anyhow!("Unable to open input file {path:?}: {error}"))?,
    ))
    .peekable();
    let mut state = State::FileStart;
    let mut current_record = FastaRecord {
        id: String::new(),
        comment: String::new(),
        sequence_handle: String::new(),
    };
    let mut records = Vec::<FastaRecord<String>>::new();

    'parser: loop {
        match state {
            State::FileStart => {
                let mut newline = true;

                'find_first_record: loop {
                    match input.next() {
                        Some(result) => match result? {
                            Character::Newline => newline = true,
                            Character::RecordStart => {
                                if newline {
                                    state = State::ParseId;
                                    break 'find_first_record;
                                } else {
                                    return Err(anyhow!(
                                        "First fasta record is not preceded by a newline character"
                                    ));
                                }
                            }
                            Character::Other(c) => {
                                newline = false;
                                if !c.is_whitespace() {
                                    return Err(anyhow!(
                                        "Found non-whitespace character before first fasta record: {c}"
                                    ));
                                }
                            }
                        },
                        None => {
                            return Err(anyhow!("Input file {path:?} contains no fasta record"));
                        }
                    }
                }
            }
            State::ParseId => 'collect_id: loop {
                match input.next() {
                    Some(result) => match result? {
                        Character::Newline => {
                            state = State::ParseSequence;
                            break 'collect_id;
                        }
                        Character::RecordStart => current_record.id.push('>'),
                        Character::Other(c) => {
                            if c.is_whitespace() {
                                state = State::ParseComment;
                                break 'collect_id;
                            } else {
                                current_record.id.push(c);
                            }
                        }
                    },
                    None => {
                        records.push(current_record);
                        break 'parser;
                    }
                }
            },
            State::ParseComment => 'collect_comment: loop {
                match input.next() {
                    Some(result) => match result? {
                        Character::Newline => {
                            state = State::ParseSequence;
                            break 'collect_comment;
                        }
                        Character::RecordStart => current_record.comment.push('>'),
                        Character::Other(c) => current_record.comment.push(c),
                    },
                    None => {
                        records.push(current_record);
                        break 'parser;
                    }
                }
            },
            State::ParseSequence => {
                let mut newline = true;

                'collect_sequence: loop {
                    match input.next() {
                        Some(result) => match result? {
                            Character::Newline => newline = true,
                            Character::RecordStart => {
                                if newline {
                                    records.push(current_record);
                                    current_record = FastaRecord {
                                        id: String::new(),
                                        comment: String::new(),
                                        sequence_handle: String::new(),
                                    };
                                    state = State::ParseId;
                                    break 'collect_sequence;
                                } else {
                                    current_record.sequence_handle.push('>');
                                    newline = false;
                                }
                            }
                            Character::Other(c) => {
                                current_record.sequence_handle.push(c);
                                newline = false;
                            }
                        },
                        None => {
                            records.push(current_record);
                            break 'parser;
                        }
                    }
                }
            }
        }
    }

    if pair {
        if records.len() != 2 {
            Err(anyhow!(
                "Expected paired fasta file with two records, but found {} records",
                records.len()
            ))
        } else {
            Ok(records)
        }
    } else if records.len() != 1 {
        Err(anyhow!(
            "Expected single-record fasta file, but found {} records",
            records.len()
        ))
    } else {
        Ok(records)
    }
}

enum Character {
    Newline,
    RecordStart,
    Other(char),
}

struct CharacterIterator<Reader: Read + ?Sized> {
    reader: BufReader<Reader>,
}

impl<Reader: Read> CharacterIterator<Reader> {
    fn new(reader: BufReader<Reader>) -> Self {
        Self { reader }
    }
}

impl<Reader: Read + ?Sized> Iterator for CharacterIterator<Reader> {
    type Item = Result<Character>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader
            .read_char_raw()
            .map(|result| {
                result.map(|c| {
                    if c == '\n' || c == '\r' {
                        Character::Newline
                    } else if c == '>' {
                        Character::RecordStart
                    } else {
                        Character::Other(c)
                    }
                })
            })
            .map_err(|error| anyhow!("Error reading character from fasta input file: {error}"))
            .transpose()
    }
}
