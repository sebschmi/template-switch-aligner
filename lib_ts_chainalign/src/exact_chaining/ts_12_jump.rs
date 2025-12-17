use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{
        Alignment, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
        ts_kind::TsDescendant,
    },
    anchors::secondary::SecondaryAnchor,
    costs::AlignmentCosts,
    exact_chaining::ts_12_jump::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct Ts12JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    alignment_costs: &'alignment_costs AlignmentCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
}

impl<'sequences, 'alignment_costs, 'rc_fn, Cost: AStarCost>
    Ts12JumpAligner<'sequences, 'alignment_costs, 'rc_fn, Cost>
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
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) -> (Cost, Alignment) {
        assert!(start.is_primary());
        assert!(end.is_secondary());

        // Enfore non-match if there is a gap between the anchors in the descendant.
        // This is to match the lower-bound computation.
        // It also discourages chains to deviate from the alignment geometry boundaries.
        let descendant_start = match end.ts_kind().unwrap().descendant {
            TsDescendant::Seq1 => start.primary_ordinate_a(),
            TsDescendant::Seq2 => start.primary_ordinate_b(),
        }
        .unwrap();
        let enforce_non_match = descendant_start != end.secondary_ordinate_descendant().unwrap();

        let context = Context::new(
            self.alignment_costs,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
            self.max_match_run,
        );
        let mut a_star = AStar::<_>::new_with_buffers(context, self.a_star_buffers.take().unwrap());

        a_star.initialise();
        let (cost, alignment) = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => (cost.0, a_star.reconstruct_path().into()),
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => (Cost::max_value(), Vec::new().into()),
        };

        a_star.search_until(|_, node| node.cost > cost);
        additional_secondary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| {
                    node.identifier.coordinates().is_secondary()
                        && (node.identifier.has_non_match() || !enforce_non_match)
                })
                .map(|node| {
                    (
                        SecondaryAnchor::new_from_start(&node.identifier.coordinates()),
                        node.cost,
                    )
                }),
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }
}
