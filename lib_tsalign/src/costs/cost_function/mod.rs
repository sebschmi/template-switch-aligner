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
