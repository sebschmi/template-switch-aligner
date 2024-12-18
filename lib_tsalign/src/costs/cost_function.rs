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

impl<SourceType: Bounded> CostFunction<SourceType> {
    /// Constructs a cost function that returns `Cost::MAX` for all input values.
    pub fn new_max() -> Self {
        Self {
            function: vec![(SourceType::min_value(), Cost::MAX)],
        }
    }
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

#[cfg(test)]
mod tests {
    use generic_a_star::cost::Cost;

    use super::CostFunction;

    #[test]
    #[expect(clippy::reversed_empty_ranges)]
    fn min() {
        let cost_function = CostFunction::try_from(vec![
            (2, Cost::from(100)),
            (3, Cost::from(1)),
            (4, Cost::from(2)),
            (6, Cost::from(1)),
            (8, Cost::from(3)),
            (70, Cost::from(2)),
            (100, Cost::from(100)),
        ])
        .unwrap();

        assert_eq!(cost_function.min(0..2), None);
        assert_eq!(cost_function.min(1..2), None);
        assert_eq!(cost_function.min(2..2), None);
        assert_eq!(cost_function.min(4..2), None);
        assert_eq!(cost_function.min(4..=2), None);
        assert_eq!(cost_function.min(3..=2), None);
        assert_eq!(cost_function.min(2..=2), Some(100.into()));
        assert_eq!(cost_function.min(3..=3), Some(1.into()));
        assert_eq!(cost_function.min(4..=4), Some(2.into()));
        assert_eq!(cost_function.min(5..=5), Some(2.into()));
        assert_eq!(cost_function.min(6..=6), Some(1.into()));
        assert_eq!(cost_function.min(2..3), Some(100.into()));
        assert_eq!(cost_function.min(3..4), Some(1.into()));
        assert_eq!(cost_function.min(4..5), Some(2.into()));
        assert_eq!(cost_function.min(5..6), Some(2.into()));
        assert_eq!(cost_function.min(6..7), Some(1.into()));
        assert_eq!(cost_function.min(22..=33), Some(3.into()));
        assert_eq!(cost_function.min(22..33), Some(3.into()));

        assert_eq!(cost_function.min(..), Some(1.into()));
        assert!([
            (0.., Some(1)),
            (1.., Some(1)),
            (2.., Some(1)),
            (3.., Some(1)),
            (4.., Some(1)),
            (5.., Some(1)),
            (6.., Some(1)),
            (7.., Some(1)),
            (8.., Some(2)),
            (9.., Some(2)),
            (69.., Some(2)),
            (70.., Some(2)),
            (71.., Some(2)),
            (72.., Some(2)),
            (99.., Some(2)),
            (100.., Some(100)),
            (101.., Some(100)),
        ]
        .into_iter()
        .map(|(range, cost)| (range, cost.map(Cost::from)))
        .all(|(range, min)| cost_function.min(range) == min));

        assert!([
            (..0, None),
            (..1, None),
            (..2, None),
            (..3, Some(100)),
            (..4, Some(1)),
            (..5, Some(1)),
            (..6, Some(1)),
            (..7, Some(1)),
            (..8, Some(1)),
            (..9, Some(1)),
            (..69, Some(1)),
            (..70, Some(1)),
            (..71, Some(1)),
            (..72, Some(1)),
            (..99, Some(1)),
            (..100, Some(1)),
            (..101, Some(1)),
        ]
        .into_iter()
        .map(|(range, cost)| (range, cost.map(Cost::from)))
        .all(|(range, min)| cost_function.min(range) == min));

        assert!([
            (..=0, None),
            (..=1, None),
            (..=2, Some(100)),
            (..=3, Some(1)),
            (..=4, Some(1)),
            (..=5, Some(1)),
            (..=6, Some(1)),
            (..=7, Some(1)),
            (..=8, Some(1)),
            (..=9, Some(1)),
            (..=69, Some(1)),
            (..=70, Some(1)),
            (..=71, Some(1)),
            (..=72, Some(1)),
            (..=99, Some(1)),
            (..=100, Some(1)),
            (..=101, Some(1)),
        ]
        .into_iter()
        .map(|(range, cost)| (range, cost.map(Cost::from)))
        .all(|(range, min)| cost_function.min(range) == min));
    }
}