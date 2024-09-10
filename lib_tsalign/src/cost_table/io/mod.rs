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
    multi::{count, many0, many1},
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

fn parse_plain<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, GapAffineAlignmentCostTable<AlphabetType>> {
    let (input, name) = opt(parse_name)(input)?;
    let (input, substitution_cost_table) = parse_substitution_cost_table::<AlphabetType>(input)?;
    let (input, gap_open_cost_vector) = parse_gap_open_cost_vector::<AlphabetType>(input)?;
    let (input, gap_extend_cost_vector) = parse_gap_extend_cost_vector::<AlphabetType>(input)?;

    let name = name.unwrap_or("").to_string();

    let cost_table = GapAffineAlignmentCostTable {
        name,
        substitution_cost_table,
        gap_open_cost_vector,
        gap_extend_cost_vector,
        phantom_data: Default::default(),
    };

    Ok((input, cost_table))
}

fn parse_name(input: &str) -> IResult<&str, &str> {
    let input = skip_any_whitespace(input)?;
    let input = char('#')(input)?.0;
    let input = many0(satisfy(is_whitespace))(input)?.0;
    take_till1(is_any_whitespace)(input)
}

fn parse_substitution_cost_table<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("SubstitutionCostTable")(input)?.0;

    // First row gives the order of the characters in the columns
    let (input, column_character_order) =
        parse_substitution_cost_table_first_row::<AlphabetType>(input)?;

    // Next is a fancy separator line
    let input = tuple((
        parse_any_whitespace,
        many1(tag("-")),
        tag("+"),
        many1(tag("-")),
    ))(input)?
    .0;

    // Then we have the rows
    let (input, mut rows) = count(
        parse_substitution_cost_table_row::<AlphabetType>,
        AlphabetType::SIZE,
    )(input)?;
    rows.sort_unstable_by_key(|(character, _)| character.index());

    // And finally we can write everything into a matrix, ensuring that everything is in the correct order
    let mut substitution_cost_table = Vec::with_capacity(AlphabetType::SIZE * AlphabetType::SIZE);
    for (_, row) in rows.iter() {
        for column in 0..AlphabetType::SIZE {
            let column_index = column_character_order
                .iter()
                .position(|character| character.index() == column)
                .unwrap();
            substitution_cost_table.push(row[column_index]);
        }
    }

    Ok((input, substitution_cost_table))
}

fn parse_substitution_cost_table_first_row<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, Vec<AlphabetType::CharacterType>> {
    let input = skip_any_whitespace(input)?;
    let input = tag("|")(input)?.0;
    let (input, characters) = count(
        preceded(
            parse_whitespace,
            parse_alphabet_character::<AlphabetType::CharacterType>,
        ),
        AlphabetType::SIZE,
    )(input)?;

    let mut sorted_characters = characters.clone();
    sorted_characters.sort();
    sorted_characters.dedup();

    if characters.len() != AlphabetType::SIZE || sorted_characters.len() != AlphabetType::SIZE {
        Err(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))
    } else {
        Ok((input, characters))
    }
}

fn parse_substitution_cost_table_row<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, (AlphabetType::CharacterType, Vec<Cost>)> {
    let input = skip_any_whitespace(input)?;
    let (input, character) = parse_alphabet_character(input)?;
    let input = skip_whitespace(input)?;
    let input = tag("|")(input)?.0;
    let (input, cost_vector) = count(
        preceded(parse_whitespace, nom::character::complete::u64),
        AlphabetType::SIZE,
    )(input)?;
    let cost_vector = cost_vector.into_iter().map(Cost::from).collect();
    Ok((input, (character, cost_vector)))
}

fn parse_gap_open_cost_vector<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("GapOpenCostVector")(input)?.0;

    parse_cost_vector::<AlphabetType>(input)
}

fn parse_gap_extend_cost_vector<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("GapExtendCostVector")(input)?.0;

    parse_cost_vector::<AlphabetType>(input)
}

fn parse_cost_vector<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    let (input, index_row) = parse_cost_vector_index_row::<AlphabetType>(input)?;
    let (input, value_row) = parse_cost_vector_value_row::<AlphabetType>(input)?;

    let mut combined_row: Vec<_> = index_row.into_iter().zip(value_row).collect();
    combined_row.sort_unstable_by_key(|(character, _)| character.index());
    let cost_vector = combined_row.into_iter().map(|(_, cost)| cost).collect();
    Ok((input, cost_vector))
}

fn parse_cost_vector_index_row<AlphabetType: Alphabet>(
    input: &str,
) -> IResult<&str, Vec<AlphabetType::CharacterType>> {
    let input = skip_any_whitespace(input)?;
    let (input, characters) = count(
        preceded(
            parse_whitespace,
            parse_alphabet_character::<AlphabetType::CharacterType>,
        ),
        AlphabetType::SIZE,
    )(input)?;

    let mut sorted_characters = characters.clone();
    sorted_characters.sort();
    sorted_characters.dedup();

    if characters.len() != AlphabetType::SIZE || sorted_characters.len() != AlphabetType::SIZE {
        Err(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))
    } else {
        Ok((input, characters))
    }
}

fn parse_cost_vector_value_row<AlphabetType: Alphabet>(input: &str) -> IResult<&str, Vec<Cost>> {
    let input = skip_any_whitespace(input)?;
    let (input, cost_vector) = count(
        preceded(parse_whitespace, nom::character::complete::u64),
        AlphabetType::SIZE,
    )(input)?;
    let cost_vector = cost_vector.into_iter().map(Cost::from).collect();
    Ok((input, cost_vector))
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
