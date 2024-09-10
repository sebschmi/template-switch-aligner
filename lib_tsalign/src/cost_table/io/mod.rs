use super::GapAffineAlignmentCostTable;
use crate::{
    cost::Cost,
    error::{Error, Result},
};

use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};
use nom::{
    bytes::complete::{tag, take, take_till1},
    character::complete::{char, satisfy},
    combinator::opt,
    multi::{many0, many1},
    sequence::{preceded, tuple},
    IResult,
};

use std::io::{Read, Write};

impl<AlphabetType: Alphabet> GapAffineAlignmentCostTable<AlphabetType> {
    pub fn read_plain(mut reader: impl Read) -> Result<Self> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;

        parse_plain(&input)
            .map_err(|error| match error {
                nom::Err::Incomplete(needed) => Error::ParserIncomplete(needed),
                nom::Err::Error(error) | nom::Err::Failure(error) => Error::Parser {
                    input: error.input.to_string(),
                    kind: error.code,
                },
            })
            .map(|(_, output)| output)
    }

    #[expect(unused)]
    pub fn write_plain(&self, writer: impl Write) -> Result<()> {
        todo!()
    }
}

#[expect(unused)]
fn parse_plain<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, GapAffineAlignmentCostTable<AlphabetType>> {
    let (input, name) = opt(parse_name)(input)?;
    let (input, substitution_table) = parse_substitution_table::<AlphabetType>(input)?;

    let name = name.unwrap_or("").to_string();

    todo!()
}

fn parse_name(input: &str) -> IResult<&str, &str> {
    let input = skip_any_whitespace(input)?;
    let input = char('#')(input)?.0;
    let input = many0(satisfy(is_whitespace))(input)?.0;
    take_till1(is_any_whitespace)(input)
}

#[expect(unused)]
fn parse_substitution_table<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("SubstitutionTable")(input)?.0;

    // First row gives the order of the characters in the columns
    let (input, column_character_order) =
        parse_substitution_table_first_row::<AlphabetType>(input)?;

    // Next is a fancy separator line
    let input = tuple((
        parse_any_whitespace,
        many1(tag("-")),
        tag("+"),
        many1(tag("-")),
    ))(input)?
    .0;

    // Then we have the rows
    todo!()
}

fn parse_substitution_table_first_row<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, Vec<AlphabetType::CharacterType>> {
    let input = skip_any_whitespace(input)?;
    let input = tag("|")(input)?.0;
    let (input, characters) = many1(preceded(
        parse_whitespace,
        parse_alphabet_character::<AlphabetType::CharacterType>,
    ))(input)?;

    if characters.len() == AlphabetType::SIZE {
        Ok((input, characters))
    } else {
        Err(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))
    }
}

fn parse_alphabet_character<CharacterType: AlphabetCharacter>(
    input: &str,
) -> IResult<&str, CharacterType> {
    let (input, character) = take(1usize)(input)?;
    let character = character.chars().next().unwrap();
    let character = CharacterType::try_from(character).map_err(|_| {
        nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        })
    })?;
    Ok((input, character))
}

fn parse_whitespace(input: &str) -> IResult<&str, ()> {
    skip_whitespace(input).map(|input| (input, ()))
}

fn parse_any_whitespace(input: &str) -> IResult<&str, ()> {
    skip_any_whitespace(input).map(|input| (input, ()))
}

fn skip_whitespace(input: &str) -> std::result::Result<&str, nom::Err<nom::error::Error<&str>>> {
    many0(satisfy(is_whitespace))(input).map(|(input, _)| input)
}

fn skip_any_whitespace(
    input: &str,
) -> std::result::Result<&str, nom::Err<nom::error::Error<&str>>> {
    many0(satisfy(is_any_whitespace))(input).map(|(input, _)| input)
}

fn is_any_whitespace(c: char) -> bool {
    c.is_whitespace() || c == '\n' || c == '\r'
}

fn is_whitespace(c: char) -> bool {
    c.is_whitespace()
}
