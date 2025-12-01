use generic_a_star::{AStar, AStarNode, cost::AStarCost};
use ndarray::{Array1, Array2};

use crate::{costs::GapAffineCosts, lower_bounds::gap_affine::algo::Context};

mod algo;
#[cfg(test)]
mod tests;

pub struct GapAffineLowerBounds<Cost> {
    lower_bounds: Array2<Cost>,
    variable_gap2_lower_bounds: Array1<Cost>,
}

impl<Cost: AStarCost> GapAffineLowerBounds<Cost> {
    pub fn new(max_n: usize, max_match_run: u32, cost_table: &GapAffineCosts<Cost>) -> Self {
        Self::compute(max_n, max_match_run, cost_table, false)
    }

    pub(super) fn new_allow_all_matches(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
    ) -> Self {
        Self::compute(max_n, max_match_run, cost_table, true)
    }

    fn compute(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineCosts<Cost>,
        allow_all_match_run: bool,
    ) -> Self {
        let mut lower_bounds = Array2::from_elem((max_n + 1, max_n + 1), Cost::max_value());
        lower_bounds[[0, 0]] = Cost::zero();
        let context = Context::new(cost_table, max_match_run, max_n);
        let mut a_star = AStar::new(context);
        a_star.initialise();
        a_star.search_until(|_, node| {
            if node.identifier.has_non_match || allow_all_match_run {
                let lower_bound = &mut lower_bounds[[node.identifier.a, node.identifier.b]];
                *lower_bound = (*lower_bound).min(node.cost().0);
            }
            false
        });
        let variable_gap2_lower_bounds = Array1::from_iter((0..=max_n).map(|gap1| {
            (0..=max_n)
                .map(|gap2| lower_bounds[[gap1, gap2]])
                .min()
                .unwrap()
        }));

        Self {
            lower_bounds,
            variable_gap2_lower_bounds,
        }
    }
}

impl<Cost: Copy> GapAffineLowerBounds<Cost> {
    /// A lower bound of the cost for chaining two anchors with the given gaps.
    /// The lower bound is symmetric, so the order of the gaps does not matter.
    pub fn lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.lower_bounds[[gap1, gap2]]
    }

    /// A lower bound of the cost for chaining two anchors with only one specified gap length.
    pub fn variable_gap2_lower_bound(&self, gap: usize) -> Cost {
        self.variable_gap2_lower_bounds[[gap]]
    }
}
