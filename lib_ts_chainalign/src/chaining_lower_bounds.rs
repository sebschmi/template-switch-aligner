//! Compute lower bounds for chaining anchors with gaps.

use generic_a_star::cost::AStarCost;
use serde::{Deserialize, Serialize};

use crate::{
    chaining_lower_bounds::{gap_affine::GapAffineLowerBounds, ts_jump::TsJumpLowerBounds},
    costs::AlignmentCosts,
};

pub mod gap_affine;
pub mod ts_jump;

#[derive(Serialize, Deserialize)]
pub struct ChainingLowerBounds<Cost> {
    primary: GapAffineLowerBounds<Cost>,
    secondary: GapAffineLowerBounds<Cost>,
    jump: TsJumpLowerBounds<Cost>,
    alignment_costs: AlignmentCosts<Cost>,
    max_match_run: u32,
}

impl<Cost: AStarCost> ChainingLowerBounds<Cost> {
    /// Compute chaining lower bounds.
    ///
    /// * `max_n` is the maximum sequence length that the lower bounds should support.
    /// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
    ///   Set this to `k-1`, if the anchors are `k`-mers.
    /// * `alignment_costs` is the cost function for the alignment.
    pub fn new(max_n: usize, max_match_run: u32, alignment_costs: AlignmentCosts<Cost>) -> Self {
        Self {
            primary: GapAffineLowerBounds::new(
                max_n,
                max_match_run,
                &alignment_costs.primary_costs,
            ),
            secondary: GapAffineLowerBounds::new(
                max_n,
                max_match_run,
                &alignment_costs.secondary_costs,
            ),
            jump: TsJumpLowerBounds::new(max_n, max_match_run, &alignment_costs),
            alignment_costs,
            max_match_run,
        }
    }
}

impl<Cost: Copy> ChainingLowerBounds<Cost> {
    pub fn primary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.primary.lower_bound(gap1, gap2)
    }

    pub fn secondary_lower_bound(&self, gap1: usize, gap2: usize) -> Cost {
        self.secondary.lower_bound(gap1, gap2)
    }

    pub fn jump_12_lower_bound(&self, descendant_gap: usize) -> Cost {
        self.jump.lower_bound_12(descendant_gap)
    }

    pub fn jump_34_lower_bound(&self, descendant_gap: usize) -> Cost {
        self.jump.lower_bound_34(descendant_gap)
    }
}

impl<Cost> ChainingLowerBounds<Cost> {
    pub fn primary(&self) -> &GapAffineLowerBounds<Cost> {
        &self.primary
    }

    pub fn secondary(&self) -> &GapAffineLowerBounds<Cost> {
        &self.secondary
    }

    pub fn jump(&self) -> &TsJumpLowerBounds<Cost> {
        &self.jump
    }

    pub fn alignment_costs(&self) -> &AlignmentCosts<Cost> {
        &self.alignment_costs
    }

    pub fn max_match_run(&self) -> u32 {
        self.max_match_run
    }
}
