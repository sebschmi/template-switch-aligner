use generic_a_star::{AStar, AStarNode, cost::AStarCost};
use ndarray::Array2;

use crate::lower_bounds::gap_affine::algo::Context;

mod algo;
#[cfg(test)]
mod tests;

pub struct GapAffineLowerBounds<Cost> {
    max_n: usize,
    lower_bounds: Array2<Cost>,
}

pub struct GapAffineLowerBoundCostTable<Cost> {
    pub substitution: Cost,
    pub gap_open: Cost,
    pub gap_extend: Cost,
}

impl<Cost: AStarCost> GapAffineLowerBounds<Cost> {
    #[expect(dead_code)]
    pub fn new(
        max_n: usize,
        max_match_run: u32,
        cost_table: &GapAffineLowerBoundCostTable<Cost>,
    ) -> Self {
        let mut lower_bounds = Array2::<Cost>::from_elem((max_n + 1, max_n + 1), Cost::max_value());
        lower_bounds[[0, 0]] = Cost::zero();
        let context = Context::new(cost_table, max_match_run, max_n);
        let mut a_star = AStar::new(context);
        a_star.initialise();
        a_star.search_until(|_, node| {
            if node.identifier.has_non_match {
                let lower_bound = &mut lower_bounds[[node.identifier.a, node.identifier.b]];
                *lower_bound = (*lower_bound).min(node.cost().0);
            }
            false
        });

        Self {
            max_n,
            lower_bounds,
        }
    }
}

impl<Cost> GapAffineLowerBounds<Cost> {
    #[expect(dead_code)]
    pub fn max_n(&self) -> usize {
        self.max_n
    }
}

impl<Cost: Copy> GapAffineLowerBounds<Cost> {
    #[expect(dead_code)]
    pub fn lower_bound(&self, a: usize, b: usize) -> Cost {
        self.lower_bounds[[a, b]]
    }
}
