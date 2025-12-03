use generic_a_star::{AStar, AStarResult, cost::AStarCost};

use crate::{
    alignment::{Alignment, coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    costs::AlignmentCosts,
    exact_chaining::ts_34_jump::algo::Context,
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts34JumpAlignment<Cost> {
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    alignment: Alignment,
    cost: Cost,
}

impl<Cost: AStarCost> Ts34JumpAlignment<Cost> {
    pub fn new(
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        sequences: &AlignmentSequences,
        cost_table: &AlignmentCosts<Cost>,
        rc_fn: &dyn Fn(u8) -> u8,
        max_match_run: u32,
    ) -> Self {
        assert!(start.is_secondary());
        assert!(end.is_primary());

        let context = Context::new(cost_table, sequences, rc_fn, start, end, max_match_run);
        let mut a_star = AStar::new(context);
        a_star.initialise();
        match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => Self {
                start,
                end,
                alignment: a_star.reconstruct_path().into(),
                // The TS base cost is applied at the 12-jump, but we anyways apply it in this algorithm to make it label-setting if the base cost is non-zero.
                // But since the 34-jump has zero cost, we subtract it again.
                cost: cost.0 - cost_table.ts_base_cost,
            },
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => Self {
                start,
                end,
                alignment: Vec::new().into(),
                cost: Cost::max_value(),
            },
        }
    }
}

impl<Cost> Ts34JumpAlignment<Cost> {
    pub fn start(&self) -> AlignmentCoordinates {
        self.start
    }

    pub fn end(&self) -> AlignmentCoordinates {
        self.end
    }

    pub fn alignment(&self) -> &Alignment {
        &self.alignment
    }
}

impl<Cost: Copy> Ts34JumpAlignment<Cost> {
    pub fn cost(&self) -> Cost {
        self.cost
    }
}
