use std::ops::{Add, Bound, RangeBounds, Sub};

use num_traits::{Bounded, One};

use crate::error::Error;

use super::cost::Cost;

pub mod io;

/// A step-wise cost funtion.
///
/// The function is represented as a list of points sorted by the input coordinate.
/// Its domain starts at the input coordinate of the first point.
/// The points mark where it steps up.
/// For example, if a cost function `f` is represented by `[(0, 1), (2, 3)]`, then
/// * `f(x)` panics for `x < 0`,
/// * `f(x) = 1` for `0 <= x < 2`, and
/// * `f(x) = 3` for `2 <= x`.
///
/// The function can be evaluated via its [`evaluate`](Self::evaluate) function.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct CostFunction<SourceType> {
    function: Vec<(SourceType, Cost)>,
}

impl<SourceType: Clone + Ord> CostFunction<SourceType> {
    /// Evaluate the cost function at position `input`.
    ///
    /// **Panics** if the given input is before the first entry in the cost function.
    pub fn evaluate(&self, input: &SourceType) -> Cost {
        match self
            .function
            .binary_search_by_key(input, |(input, _)| input.clone())
        {
            Ok(index) => self.function[index].1,
            Err(index) => self.function[index - 1].1,
        }
    }

    pub fn minimum_finite_input(&self) -> Option<SourceType> {
        self.function
            .iter()
            .filter_map(|(input, cost)| {
                if *cost < Cost::MAX {
                    Some(input.clone())
                } else {
                    None
                }
            })
            .next()
    }
}

impl<
        SourceType: Clone + Ord + Bounded + One + Add<Output = SourceType> + Sub<Output = SourceType>,
    > CostFunction<SourceType>
{
    pub fn min(&self, range: impl RangeBounds<SourceType>) -> Option<Cost> {
        let is_not_empty = match (range.start_bound(), range.end_bound()) {
            (Bound::Included(start), Bound::Included(end)) => start <= end,
            (Bound::Included(start), Bound::Excluded(end)) => start < end,
            (Bound::Excluded(start), Bound::Included(end)) => start < end,
            (Bound::Excluded(start), Bound::Excluded(end)) => {
                start.clone() + SourceType::one() < end.clone()
            }
            (Bound::Excluded(start), Bound::Unbounded) => start != &SourceType::max_value(),
            (Bound::Unbounded, Bound::Excluded(end)) => end != &SourceType::min_value(),
            (Bound::Included(_), Bound::Unbounded) => true,
            (Bound::Unbounded, Bound::Included(_)) => true,
            (Bound::Unbounded, Bound::Unbounded) => true,
        };

        if !is_not_empty {
            return None;
        }

        self.function
            .windows(2)
            .filter_map(|window| {
                let (first_input, first_cost) = &window[0];
                let (next_input, _) = &window[1];
                let last_input = next_input.clone() - SourceType::one();

                let first_left_of_end = match range.end_bound() {
                    Bound::Included(end) => first_input <= end,
                    Bound::Excluded(end) => first_input < end,
                    Bound::Unbounded => true,
                };
                let last_right_of_start = match range.start_bound() {
                    Bound::Included(start) => start <= &last_input,
                    Bound::Excluded(start) => start < &last_input,
                    Bound::Unbounded => true,
                };

                if first_left_of_end && last_right_of_start {
                    Some(first_cost)
                } else {
                    None
                }
            })
            .chain(
                self.function
                    .last()
                    .iter()
                    .filter_map(|(first_input, cost)| {
                        let first_left_of_end = match range.end_bound() {
                            Bound::Included(end) => first_input <= end,
                            Bound::Excluded(end) => first_input < end,
                            Bound::Unbounded => true,
                        };

                        if first_left_of_end {
                            Some(cost)
                        } else {
                            None
                        }
                    }),
            )
            .min()
            .copied()
    }
}

impl<SourceType: Clone + Ord + From<u8> + Sub<Output = SourceType>> CostFunction<SourceType> {
    pub fn maximum_finite_input(&self) -> Option<SourceType> {
        let last_finite_index = self
            .function
            .iter()
            .enumerate()
            .rev()
            .find_map(|(index, (_, cost))| if *cost < Cost::MAX { Some(index) } else { None })?;
        let infinite_index = last_finite_index + 1;

        if infinite_index == self.function.len() {
            return None;
        }

        Some(self.function[infinite_index].0.clone() - 1.into())
    }
}

impl<SourceType: Ord> TryFrom<Vec<(SourceType, Cost)>> for CostFunction<SourceType> {
    type Error = Error;
    fn try_from(function: Vec<(SourceType, Cost)>) -> Result<Self, Self::Error> {
        for (index, window) in function.windows(2).enumerate() {
            if window[0].0 >= window[1].0 {
                return Err(Error::CostFunctionIndexNotIncreasing { index: index + 1 });
            }
        }

        Ok(Self { function })
    }
}

impl<SourceType> From<CostFunction<SourceType>> for Vec<(SourceType, Cost)> {
    fn from(value: CostFunction<SourceType>) -> Self {
        value.function
    }
}
