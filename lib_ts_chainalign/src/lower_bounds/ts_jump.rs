use generic_a_star::cost::AStarCost;
use ndarray::{Array1, Array3};

use crate::{costs::AlignmentCosts, lower_bounds::gap_affine::GapAffineLowerBounds};

pub struct TsJumpLowerBounds<Cost> {
    primary_lower_bounds: GapAffineLowerBounds<Cost>,
    secondary_lower_bounds: GapAffineLowerBounds<Cost>,
    lower_bounds_12: Array1<Cost>,
    lower_bounds_1234: Array3<Cost>,
}

impl<Cost: AStarCost> TsJumpLowerBounds<Cost> {
    pub fn new(max_n: usize, max_match_run: u32, cost_table: &AlignmentCosts<Cost>) -> Self {
        let primary_lower_bounds =
            GapAffineLowerBounds::new(max_n, max_match_run, &cost_table.primary_costs);
        let secondary_lower_bounds =
            GapAffineLowerBounds::new(max_n, max_match_run, &cost_table.secondary_costs);

        // This way of calculating the lower bound for the 12-jump does not take the shape limits of the template switch into account.
        // However, most of the time these limits are gonna be big, so they should not have a big impact on the lower bound.
        let mut lower_bounds_12 = Array1::<Cost>::from_elem(max_n + 1, Cost::max_value());
        for primary_descendant_gap in 0..=max_n {
            for secondary_descendant_gap in 0..=max_n - primary_descendant_gap {
                let lower_bound = primary_lower_bounds
                    .variable_gap2_lower_bound(primary_descendant_gap)
                    + cost_table.ts_base_cost
                    + secondary_lower_bounds.variable_gap2_lower_bound(secondary_descendant_gap);
                lower_bounds_12[[primary_descendant_gap + secondary_descendant_gap]] =
                    lower_bounds_12[[primary_descendant_gap + secondary_descendant_gap]]
                        .min(lower_bound);
            }
        }

        todo!()
    }
}

impl<Cost: Copy> TsJumpLowerBounds<Cost> {
    /// A lower bound of the cost for chaining two primary anchors with the given gaps.
    /// The lower bound is symmetric, so the order of the gaps does not matter.
    #[expect(dead_code)]
    pub fn primary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.primary_lower_bounds.lower_bound(gap1, gap2)
    }

    /// A lower bound of the cost for chaining two secondary anchors with the given gaps.
    /// The lower bound is symmetric, so the order of the gaps does not matter.
    #[expect(dead_code)]
    pub fn secondary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.secondary_lower_bounds.lower_bound(gap1, gap2)
    }

    /// A lower bound of the cost for chaining a primary anchor with a secondary anchor.
    /// As the ancestor gap is determined by the 34-jump which is not know when the 12-jump is evaluated,
    /// this lower bound only depends on the descendant gap.
    ///
    /// This lower bound takes the template switch base cost into account.
    #[expect(dead_code)]
    pub fn lower_bound_12(&self, descendant_gap: usize) -> Cost {
        self.lower_bounds_12[[descendant_gap]]
    }

    /// A lower bound of the cost of the jump chainings of a template switch.
    /// This is a bound to chaining a primary anchor with a secondary anchor for the 12-jump, then chaining some secondary anchors (possibly none),
    /// and finally chaining the last secondary anchor with a primary anchor for the 34-jump.
    ///
    /// **Note:** This lower bound supersedes the 12-jump lower bound, so adding both together would be incorrect.
    #[expect(dead_code)]
    pub fn lower_bound_1234(
        &self,
        ancestor_gap: usize,
        descendant_gap1: usize,
        descendant_gap2: usize,
    ) -> Cost {
        self.lower_bounds_1234[[ancestor_gap, descendant_gap1, descendant_gap2]]
    }
}
