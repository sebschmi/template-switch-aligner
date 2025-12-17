use generic_a_star::{AStar, AStarBuffers, AStarResult, cost::AStarCost};

use crate::{
    alignment::{Alignment, coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    anchors::{primary::PrimaryAnchor, secondary::SecondaryAnchor},
    costs::GapAffineCosts,
    exact_chaining::gap_affine::algo::{Context, Node},
};

mod algo;
#[cfg(test)]
mod tests;

pub struct GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost> {
    a_star_buffers: Option<AStarBuffers<Node<Cost>>>,
    sequences: &'sequences AlignmentSequences,
    cost_table: &'cost_table GapAffineCosts<Cost>,
    rc_fn: &'rc_fn dyn Fn(u8) -> u8,
    max_match_run: u32,
}

impl<'sequences, 'cost_table, 'rc_fn, Cost: AStarCost>
    GapAffineAligner<'sequences, 'cost_table, 'rc_fn, Cost>
{
    pub fn new(
        sequences: &'sequences AlignmentSequences,
        cost_table: &'cost_table GapAffineCosts<Cost>,
        rc_fn: &'rc_fn dyn Fn(u8) -> u8,
        max_match_run: u32,
    ) -> Self {
        Self {
            a_star_buffers: Some(Default::default()),
            sequences,
            cost_table,
            rc_fn,
            max_match_run,
        }
    }

    pub fn align(
        &mut self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) -> (Cost, Alignment) {
        assert!(
            start.is_primary() && end.is_primary() || start.is_secondary() && end.is_secondary()
        );

        if start == end {
            return (Cost::zero(), Vec::new().into());
        }

        let context = Context::new(
            self.cost_table,
            self.sequences,
            self.rc_fn,
            start,
            end,
            true,
            self.max_match_run,
        );
        let mut a_star = AStar::new_with_buffers(context, self.a_star_buffers.take().unwrap());

        a_star.initialise();
        let (cost, alignment) = match a_star.search() {
            AStarResult::FoundTarget { cost, .. } => {
                let alignment = a_star.reconstruct_path().into();
                (cost.0, alignment)
            }
            AStarResult::ExceededCostLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::ExceededMemoryLimit { .. } => unreachable!("Cost limit is None"),
            AStarResult::NoTarget => (Cost::max_value(), Vec::new().into()),
        };

        a_star.search_until(|_, node| node.cost > cost);
        additional_primary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| {
                    node.identifier.coordinates.is_primary() && node.identifier.has_non_match
                })
                .map(|node| {
                    (
                        PrimaryAnchor::new_from_start(&node.identifier.coordinates),
                        node.cost,
                    )
                }),
        );
        additional_secondary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| {
                    node.identifier.coordinates.is_secondary() && node.identifier.has_non_match
                })
                .map(|node| {
                    (
                        SecondaryAnchor::new_from_start(&node.identifier.coordinates),
                        node.cost,
                    )
                }),
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }
}
