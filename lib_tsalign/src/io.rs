use nom::{
    IResult, Parser,
    bytes::complete::take_till1,
    character::complete::{char, satisfy},
    multi::many0,
};

use crate::error::Error;

pub fn parse_title(input: &str) -> IResult<&str, &str> {
    let input = skip_any_whitespace(input)?;
    let input = char('#')(input)?.0;
    let input = many0(satisfy(is_whitespace)).parse(input)?.0;
    let (input, result) = take_till1(is_any_line_break)(input)?;
    Ok((input, result.trim()))
}

pub fn parse_whitespace(input: &str) -> IResult<&str, ()> {
    skip_whitespace(input).map(|input| (input, ()))
}

pub fn parse_any_whitespace(input: &str) -> IResult<&str, ()> {
    skip_any_whitespace(input).map(|input| (input, ()))
}

pub fn skip_whitespace(
    input: &str,
) -> std::result::Result<&str, nom::Err<nom::error::Error<&str>>> {
    many0(satisfy(is_whitespace))
        .parse(input)
        .map(|(input, _)| input)
}

pub fn skip_any_whitespace(
    input: &str,
) -> std::result::Result<&str, nom::Err<nom::error::Error<&str>>> {
    many0(satisfy(is_any_whitespace))
        .parse(input)
        .map(|(input, _)| input)
}

pub fn is_any_whitespace(c: char) -> bool {
    is_whitespace(c) || is_any_line_break(c)
}

pub fn is_whitespace(c: char) -> bool {
    c.is_whitespace() && !is_any_line_break(c)
}

pub fn is_any_line_break(c: char) -> bool {
    c == '\n' || c == '\r'
}

pub fn translate_nom_error(error: nom::Err<nom::error::Error<&str>>) -> Error {
    match error {
        nom::Err::Incomplete(needed) => Error::ParserIncomplete(needed),
        nom::Err::Error(error) | nom::Err::Failure(error) => Error::Parser {
            input: error.input.to_string(),
            kind: error.code,
        },
    }
}
