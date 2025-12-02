use generic_a_star::{AStar, AStarResult, cost::AStarCost};

use crate::{
    alignment::{Alignment, coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    costs::AlignmentCosts,
    exact_chaining::ts_12_jump::algo::Context,
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts12JumpAlignment<Cost> {
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
    alignment: Alignment,
    cost: Cost,
}

impl<Cost: AStarCost> Ts12JumpAlignment<Cost> {
    pub fn new(
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        sequences: &AlignmentSequences,
        cost_table: &AlignmentCosts<Cost>,
        max_match_run: u32,
    ) -> Self {
        assert!(start.is_primary());
        assert!(end.is_secondary());

        let context = Context::new(cost_table, sequences, start, end, max_match_run);
        let mut a_star = AStar::new(context);
        a_star.initialise();
        match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => Self {
                start,
                end,
                alignment: a_star.reconstruct_path().into(),
                cost: cost.0,
            },
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => {
                panic!("No TS 12-jump alignment found between the given coordinates")
            }
        }
    }
}

impl<Cost> Ts12JumpAlignment<Cost> {
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

impl<Cost: Copy> Ts12JumpAlignment<Cost> {
    pub fn cost(&self) -> Cost {
        self.cost
    }
}
