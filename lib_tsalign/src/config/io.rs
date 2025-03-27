use std::str::FromStr;

use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::AStarCost;
use log::trace;
use nom::{
    IResult,
    bytes::complete::{tag, take_while1},
    character::complete::line_ending,
    sequence::tuple,
};
use num_traits::{Bounded, PrimInt};

use crate::{
    config::BaseCost,
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    io::{parse_any_whitespace, parse_whitespace, skip_any_whitespace, translate_nom_error},
};

use super::TemplateSwitchConfig;

impl<AlphabetType: Alphabet, Cost: AStarCost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn read_plain(mut reader: impl std::io::Read) -> crate::error::Result<Self> {
        let mut input = String::new();
        reader.read_to_string(&mut input)?;
        let input = input.as_str();
        Self::parse_plain(input)
            .map(|(_, result)| result)
            .map_err(translate_nom_error)
    }

    fn parse_plain(input: &str) -> IResult<&str, Self> {
        trace!("Parsing limits");
        let (input, ()) = parse_specific_name("Limits")(input)?;
        let (input, left_flank_length) = parse_specific_equals_value("left_flank_length")(input)?;
        let (input, right_flank_length) = parse_specific_equals_value("right_flank_length")(input)?;

        trace!("Parsing base costs");
        let (input, ()) = parse_specific_name("Base Cost")(input)?;
        let (input, rr_cost) = parse_specific_equals_value("rr_cost")(input)?;
        let (input, rq_cost) = parse_specific_equals_value("rq_cost")(input)?;
        let (input, qr_cost) = parse_specific_equals_value("qr_cost")(input)?;
        let (input, qq_cost) = parse_specific_equals_value("qq_cost")(input)?;

        trace!("Parsing jump costs");
        let (input, ()) = parse_specific_name("Jump Costs")(input)?;
        let (input, offset_costs) = parse_named_cost_function("Offset")(input)?;
        let (input, length_costs) = parse_named_cost_function("Length")(input)?;
        let (input, length_difference_costs) =
            parse_named_cost_function("LengthDifference")(input)?;

        trace!("Parsing primary edit costs");
        let (input, primary_edit_costs) = parse_named_cost_table("Primary Edit Costs")(input)?;
        trace!("Parsing secondary edit costs");
        let (input, secondary_edit_costs) = parse_named_cost_table("Secondary Edit Costs")(input)?;
        trace!("Parsing left flank edit costs");
        let (input, left_flank_edit_costs) =
            parse_named_cost_table("Left Flank Edit Costs")(input)?;
        trace!("Parsing right flank edit costs");
        let (input, right_flank_edit_costs) =
            parse_named_cost_table("Right Flank Edit Costs")(input)?;

        Ok((
            input,
            Self {
                left_flank_length,
                right_flank_length,
                min_length: length_costs.minimum_finite_input().unwrap_or(usize::MAX),

                base_cost: BaseCost {
                    rr: rr_cost,
                    rq: rq_cost,
                    qr: qr_cost,
                    qq: qq_cost,
                },

                primary_edit_costs,
                secondary_edit_costs,
                left_flank_edit_costs,
                right_flank_edit_costs,

                offset_costs,
                length_costs,
                length_difference_costs,
            },
        ))
    }
}

fn parse_specific_name(name: &str) -> impl '_ + FnMut(&str) -> IResult<&str, ()> {
    move |input| {
        tuple((
            parse_any_whitespace,
            tag("#"),
            parse_whitespace,
            tag(name),
            parse_whitespace,
            line_ending,
            parse_any_whitespace,
        ))(input)
        .map(|(input, _)| (input, ()))
    }
}

fn parse_specific_equals_value<Value: FromStr>(
    identifier: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, Value> {
    move |input| {
        let (input, (actual_identifier, value)) = parse_equals_value(input)?;
        if actual_identifier == identifier {
            Ok((input, value))
        } else {
            Err(nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            }))
        }
    }
}

fn parse_equals_value<Value: FromStr>(input: &str) -> IResult<&str, (&str, Value)> {
    let input = skip_any_whitespace(input)?;
    let (input, identifier) = take_while1(|c: char| c.is_alphanumeric() || c == '_')(input)?;
    let (input, _) = tuple((parse_whitespace, tag("="), parse_whitespace))(input)?;
    let (input, value) = take_while1(|c: char| !c.is_whitespace())(input)?;

    let value = Value::from_str(value).map_err(|_| {
        nom::Err::Failure(nom::error::Error {
            input,
            code: nom::error::ErrorKind::Verify,
        })
    })?;

    Ok((input, (identifier, value)))
}

fn parse_named_cost_function<SourceType: FromStr + PrimInt, Cost: FromStr + Bounded>(
    name: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, CostFunction<SourceType, Cost>> {
    move |input| {
        let input = skip_any_whitespace(input)?;
        let (input, _) = tag(name)(input)?;
        CostFunction::parse_plain(input)
    }
}

fn parse_named_cost_table<AlphabetType: Alphabet, Cost: AStarCost>(
    name: &str,
) -> impl '_ + FnMut(&str) -> IResult<&str, GapAffineAlignmentCostTable<AlphabetType, Cost>> {
    move |input| {
        let (input, result) = GapAffineAlignmentCostTable::parse_plain(input)?;
        if result.name() == name {
            Ok((input, result))
        } else {
            Err(nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            }))
        }
    }
}
