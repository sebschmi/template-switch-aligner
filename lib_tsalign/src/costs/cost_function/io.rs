use std::{fmt::Display, io::Write, str::FromStr};

use nom::IResult;
use num_traits::{Bounded, PrimInt};

use super::CostFunction;

use crate::{
    config::io::parse_inf_value,
    error::Result,
    io::{skip_any_whitespace, skip_whitespace},
};

impl<SourceType: PrimInt + Display, Cost: Bounded + Display + Eq> CostFunction<SourceType, Cost> {
    pub fn write_plain(&self, mut writer: impl Write) -> Result<()> {
        let column_widths: Vec<_> = self
            .function
            .iter()
            .map(|(index, cost)| {
                if index == &SourceType::max_value() {
                    3
                } else if index == &SourceType::min_value() {
                    4
                } else {
                    format!("{index}").len()
                }
                .max(if cost == &Cost::max_value() {
                    3
                } else {
                    format!("{cost}").len()
                })
            })
            .collect();

        let mut once = false;
        for (&column_width, (index, _)) in column_widths.iter().zip(self.function.iter()) {
            if once {
                write!(writer, " ")?;
            } else {
                once = true;
            }

            if index == &SourceType::max_value() {
                for _ in 3..column_width {
                    write!(writer, " ")?;
                }
                write!(writer, "inf")?;
            } else if index == &SourceType::min_value() {
                for _ in 4..column_width {
                    write!(writer, " ")?;
                }
                write!(writer, "-inf")?;
            } else {
                write!(writer, "{index: >column_width$}")?;
            }
        }
        writeln!(writer)?;

        let mut once = false;
        for (&column_width, (_, cost)) in column_widths.iter().zip(self.function.iter()) {
            if once {
                write!(writer, " ")?;
            } else {
                once = true;
            }

            if cost == &Cost::max_value() {
                for _ in 3..column_width {
                    write!(writer, " ")?;
                }
                write!(writer, "inf")?;
            } else {
                write!(writer, "{cost: >column_width$}")?;
            }
        }

        Ok(())
    }
}

impl<SourceType: FromStr + PrimInt, Cost: FromStr + Bounded> CostFunction<SourceType, Cost> {
    pub(crate) fn parse_plain(input: &str) -> IResult<&str, Self> {
        let mut input = skip_any_whitespace(input)?;

        let mut indexes = Vec::new();
        while !input.starts_with(['\n', '\r']) {
            let index: SourceType;
            (input, index) = parse_inf_value(input)?;
            indexes.push(index);
            input = skip_whitespace(input)?;
        }

        input = skip_any_whitespace(input)?;

        let mut costs = Vec::new();
        while !input.starts_with(['\n', '\r']) && !input.is_empty() {
            let cost: Cost;
            (input, cost) = parse_inf_value(input)?;
            costs.push(cost);
            input = skip_whitespace(input)?;
        }

        if indexes.len() != costs.len()
            || indexes[0] != SourceType::min_value()
            || indexes.windows(2).any(|window| window[0] >= window[1])
        {
            return Err(nom::Err::Failure(nom::error::Error {
                input,
                code: nom::error::ErrorKind::Verify,
            }));
        }

        Ok((
            input,
            Self {
                function: indexes.into_iter().zip(costs).collect(),
            },
        ))
    }
}

#[cfg(test)]
mod tests {
    use generic_a_star::cost::U64Cost;

    use crate::costs::cost_function::CostFunction;

    #[test]
    fn simple_example() {
        let input = "-inf -12345 -4 -1 0 1 +2 123456 inf\n   1      2  3  4 5 6  7      8   9";
        let expected_output =
            "-inf -12345 -4 -1 0 1 2 123456 inf\n   1      2  3  4 5 6 7      8   9";
        let expected_parsing_result = CostFunction::<isize, U64Cost> {
            function: vec![
                (isize::MIN, 1u64.into()),
                (-12345, 2u64.into()),
                (-4, 3u64.into()),
                (-1, 4u64.into()),
                (0, 5u64.into()),
                (1, 6u64.into()),
                (2, 7u64.into()),
                (123456, 8u64.into()),
                (isize::MAX, 9u64.into()),
            ],
        };

        let (remaining_input, actual_parsing_result) =
            CostFunction::<isize, U64Cost>::parse_plain(input).unwrap();
        assert!(remaining_input.is_empty());

        let mut writer = Vec::new();
        actual_parsing_result.write_plain(&mut writer).unwrap();
        let output = String::from_utf8(writer).unwrap();

        assert_eq!(expected_parsing_result, actual_parsing_result);
        assert_eq!(expected_output, output);
    }
}
