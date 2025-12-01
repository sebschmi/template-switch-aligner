use generic_a_star::cost::AStarCost;

use crate::{
    costs::AlignmentCosts,
    lower_bounds::{gap_affine::GapAffineLowerBounds, ts_jump::TsJumpLowerBounds},
};

mod costs;
mod lower_bounds;

#[expect(dead_code)]
fn compute_lower_bounds<Cost: AStarCost>(
    max_n: usize,
    max_match_run: u32,
    costs: &AlignmentCosts<Cost>,
) {
    let gap_affine_lower_bounds =
        GapAffineLowerBounds::new(max_n, max_match_run, &costs.primary_costs);
    let ts_jump_lower_bounds = TsJumpLowerBounds::new(max_n, max_match_run, costs);

    // Remove dead code warnings
    gap_affine_lower_bounds.lower_bound(0, 0);
    ts_jump_lower_bounds.lower_bound_12(0);
    ts_jump_lower_bounds.lower_bound_34(0);
}
