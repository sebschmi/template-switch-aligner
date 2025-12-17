use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{Alignment, coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    costs::AlignmentCosts,
    exact_chaining::ts_34_jump::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost>
    Ts34JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
    ) -> Self {
        Self {
            a_star_buffers: Some(Default::default()),
            sequences,
            alignment_costs,
            rc_fn,
            max_match_run,
        }
    }

    pub fn align(
        &mut self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> (Cost, Alignment) {
        assert!(start.is_secondary());
        assert!(end.is_primary());

        // Do not enforce non-match, because TS geometry may cause match-only jump alignments.

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            self.max_match_run,
        );
        let mut a_star = AStar::<_>::new_with_buffers(context, self.a_star_buffers.take().unwrap());

        a_star.initialise();
        let result = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => {
                // The TS base cost is applied at the 12-jump, but we anyways apply it in this algorithm to make it label-setting if the base cost is non-zero.
                // But since the 34-jump has zero cost, we subtract it again.
                let cost = cost.0 - self.alignment_costs.ts_base_cost;
                let alignment = a_star.reconstruct_path().into();

                (cost, alignment)
            }
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => (Cost::max_value(), Vec::new().into()),
        };

        self.a_star_buffers = Some(a_star.into_buffers());

        result
    }
}
