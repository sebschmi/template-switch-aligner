use super::GapAffineAlignmentCostTable;
use crate::{
    config::io::parse_inf_value,
    error::{Error, Result},
    io::{
        parse_any_whitespace, parse_title, parse_whitespace, skip_any_whitespace, skip_whitespace,
        translate_nom_error,
    },
};

use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};
use generic_a_star::cost::AStarCost;
use log::trace;
use nom::{
    IResult, Parser,
    bytes::complete::{tag, take},
    combinator::opt,
    multi::{count, many1},
    sequence::preceded,
};

use std::{
    collections::HashMap,
    io::{Read, Write},
};

#[cfg(test)]
mod tests;

impl<AlphabetType: Alphabet, Cost: AStarCost> GapAffineAlignmentCostTable<AlphabetType, Cost> {
    pub fn read_plain_multi(mut reader: impl Read) -> Result<HashMap<String, Self>> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;
        let mut input = input.as_str();
        let mut result = HashMap::new();

        loop {
            let (next_input, table) = Self::parse_plain(input).map_err(translate_nom_error)?;
            input = next_input;
            if let Some(previous_table) = result.insert(table.name().to_string(), table) {
                return Err(Error::DuplicateCostTableName(previous_table.name));
            }

            input = skip_any_whitespace(input).map_err(translate_nom_error)?;
            if input.is_empty() {
                break;
            }
        }

        Ok(result)
    }

    pub fn read_plain(mut reader: impl Read) -> Result<Self> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;

        Self::parse_plain(&input)
            .map_err(translate_nom_error)
            .map(|(_, output)| output)
    }

    pub fn write_plain(&self, mut writer: impl Write) -> Result<()> {
        writeln!(writer, "# {}", self.name)?;
        writeln!(writer)?;

        writeln!(writer, "SubstitutionCostTable")?;

        let column_width = self
            .substitution_cost_table
            .iter()
            .map(|substitution_cost| format!("{substitution_cost}").len())
            .max()
            .unwrap();
        println!("{column_width:?}");

        write!(writer, "  |")?;
        for column_index in 0..AlphabetType::SIZE {
            let character = AlphabetType::CharacterType::from_index(column_index).unwrap();
            for _ in 0..column_width {
                write!(writer, " ")?;
            }
            write!(writer, "{character}")?;
        }
        writeln!(writer)?;

        write!(writer, "--+")?;
        for _ in 0..(usize::from(AlphabetType::SIZE) * (column_width + 1)) {
            write!(writer, "-")?;
        }
        writeln!(writer)?;

        for row_index in 0..AlphabetType::SIZE {
            let character = AlphabetType::CharacterType::from_index(row_index).unwrap();
            write!(writer, "{character} |")?;
            for column_index in 0..AlphabetType::SIZE {
                let cost = self.substitution_cost_table[usize::from(row_index)
                    * usize::from(AlphabetType::SIZE)
                    + usize::from(column_index)];
                write!(writer, " {cost: >column_width$}")?;
            }
            writeln!(writer)?;
        }
        writeln!(writer)?;

        writeln!(writer, "GapOpenCostVector")?;

        let column_width = self
            .gap_open_cost_vector
            .iter()
            .map(|gap_open_cost| format!("{gap_open_cost}").len())
            .max()
            .unwrap();

        for column_index in 0..AlphabetType::SIZE {
            let character = AlphabetType::CharacterType::from_index(column_index).unwrap();
            for _ in 0..column_width {
                write!(writer, " ")?;
            }
            write!(writer, "{character}")?;
        }
        writeln!(writer)?;

        for column_index in 0..AlphabetType::SIZE {
            let cost = self.gap_open_cost_vector[usize::from(column_index)];
            write!(writer, " {cost: >column_width$}")?;
        }
        writeln!(writer)?;
        writeln!(writer)?;

        writeln!(writer, "GapExtendCostVector")?;

        let column_width = self
            .gap_extend_cost_vector
            .iter()
            .map(|gap_extend_cost| format!("{gap_extend_cost}").len())
            .max()
            .unwrap();

        for column_index in 0..AlphabetType::SIZE {
            let character = AlphabetType::CharacterType::from_index(column_index).unwrap();
            for _ in 0..column_width {
                write!(writer, " ")?;
            }
            write!(writer, "{character}")?;
        }
        writeln!(writer)?;

        for column_index in 0..AlphabetType::SIZE {
            let cost = self.gap_extend_cost_vector[usize::from(column_index)];
            write!(writer, " {cost: >column_width$}")?;
        }
        writeln!(writer)?;

        Ok(())
    }

    pub(crate) fn parse_plain(input: &str) -> IResult<&str, Self> {
        let (input, name) = opt(parse_title).parse(input)?;
        let (input, substitution_cost_table) =
            parse_substitution_cost_table::<AlphabetType, _>(input)?;
        let (input, gap_open_cost_vector) = parse_gap_open_cost_vector::<AlphabetType, _>(input)?;
        let (input, gap_extend_cost_vector) =
            parse_gap_extend_cost_vector::<AlphabetType, _>(input)?;

        let name = name.unwrap_or("").to_string();

        let cost_table = Self {
            name,
            substitution_cost_table,
            gap_open_cost_vector,
            gap_extend_cost_vector,
            phantom_data: Default::default(),
        };

        Ok((input, cost_table))
    }
}

fn parse_substitution_cost_table<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("SubstitutionCostTable")(input)?.0;

    // First row gives the order of the characters in the columns
    let (input, column_character_order) =
        parse_substitution_cost_table_first_row::<AlphabetType>(input)?;

    // Next is a fancy separator line
    let input = (
        parse_any_whitespace,
        many1(tag("-")),
        tag("+"),
        many1(tag("-")),
    )
        .parse(input)?
        .0;

    // Then we have the rows
    let (input, mut rows) = count(
        parse_substitution_cost_table_row::<AlphabetType, _>,
        AlphabetType::SIZE.into(),
    )
    .parse(input)?;
    rows.sort_unstable_by_key(|(character, _)| character.index());

    // And finally we can write everything into a matrix, ensuring that everything is in the correct order
    let mut substitution_cost_table =
        Vec::with_capacity(usize::from(AlphabetType::SIZE) * usize::from(AlphabetType::SIZE));
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
    trace!(
        "Parsing the first row of a substitution cost table with alphabet size {}",
        AlphabetType::SIZE
    );

    let input = skip_any_whitespace(input)?;
    let input = tag("|")(input)?.0;
    let (input, characters) = count(
        preceded(
            parse_whitespace,
            parse_alphabet_character::<AlphabetType::CharacterType>,
        ),
        AlphabetType::SIZE.into(),
    )
    .parse(input)?;

    let mut sorted_characters = characters.clone();
    sorted_characters.sort();
    sorted_characters.dedup();

    if characters.len() != AlphabetType::SIZE.into()
        || sorted_characters.len() != AlphabetType::SIZE.into()
    {
        Err(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))
    } else {
        Ok((input, characters))
    }
}

fn parse_substitution_cost_table_row<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, (AlphabetType::CharacterType, Vec<Cost>)> {
    let input = skip_any_whitespace(input)?;
    let (input, character) = parse_alphabet_character(input)?;
    let input = skip_whitespace(input)?;
    let input = tag("|")(input)?.0;
    let (input, cost_vector) = count(
        preceded(parse_whitespace, parse_inf_value),
        AlphabetType::SIZE.into(),
    )
    .parse(input)?;
    let cost_vector = cost_vector.into_iter().collect();
    Ok((input, (character, cost_vector)))
}

fn parse_gap_open_cost_vector<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("GapOpenCostVector")(input)?.0;

    parse_cost_vector::<AlphabetType, _>(input)
}

fn parse_gap_extend_cost_vector<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, Vec<Cost>> {
    // Identifier
    let input = skip_any_whitespace(input)?;
    let input = tag("GapExtendCostVector")(input)?.0;

    parse_cost_vector::<AlphabetType, _>(input)
}

fn parse_cost_vector<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, Vec<Cost>> {
    let (input, index_row) = parse_cost_vector_index_row::<AlphabetType>(input)?;
    let (input, value_row) = parse_cost_vector_value_row::<AlphabetType, _>(input)?;

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
        AlphabetType::SIZE.into(),
    )
    .parse(input)?;

    let mut sorted_characters = characters.clone();
    sorted_characters.sort();
    sorted_characters.dedup();

    if characters.len() != AlphabetType::SIZE.into()
        || sorted_characters.len() != AlphabetType::SIZE.into()
    {
        Err(nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        }))
    } else {
        Ok((input, characters))
    }
}

fn parse_cost_vector_value_row<AlphabetType: Alphabet, Cost: AStarCost>(
    input: &str,
) -> IResult<&str, Vec<Cost>> {
    let input = skip_any_whitespace(input)?;
    let (input, cost_vector) = count(
        preceded(parse_whitespace, parse_inf_value),
        AlphabetType::SIZE.into(),
    )
    .parse(input)?;
    let cost_vector = cost_vector.into_iter().collect();
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
