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

        let minimum_primary_sequence_length =
            (self.sequences.primary_end().primary_ordinate_a().unwrap()
                - self.sequences.primary_start().primary_ordinate_a().unwrap())
            .min(
                self.sequences.primary_end().primary_ordinate_b().unwrap()
                    - self.sequences.primary_start().primary_ordinate_b().unwrap(),
            );
        let allow_direct_chaining =
            start == self.sequences.primary_start() || end == self.sequences.primary_end();
        let allow_all_matches = start == self.sequences.primary_start()
            && end == self.sequences.primary_end()
            && u32::try_from(minimum_primary_sequence_length).unwrap() <= self.max_match_run;

        let enforce_non_match = true;

        let context = Context::new(
            self.cost_table,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
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
        a_star.search_until_with_target_policy(|_, node| node.cost > cost, true);

        Self::fill_additional_targets(
            &a_star,
            start,
            enforce_non_match,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output lists.
    pub fn align_until_cost_limit(
        &mut self,
        start: AlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) {
        let enforce_non_match = true;
        let end = self.sequences.end(start.ts_kind());
        debug_assert!(
            start.is_primary() && end.is_primary() || start.is_secondary() && end.is_secondary()
        );

        let context = Context::new(
            self.cost_table,
            self.sequences,
            self.rc_fn,
            start,
            end,
            enforce_non_match,
            self.max_match_run,
        );
        let mut a_star = AStar::new_with_buffers(context, self.a_star_buffers.take().unwrap());
        a_star.initialise();
        a_star.search_until_with_target_policy(|_, node| node.cost > cost_limit, true);

        Self::fill_additional_targets(
            &a_star,
            start,
            enforce_non_match,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());
    }

    fn fill_additional_targets(
        a_star: &AStar<Context<Cost>>,
        start: AlignmentCoordinates,
        enforce_non_match: bool,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) {
        additional_primary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| {
                    node.identifier.coordinates.is_primary()
                        && ((node.identifier.has_non_match
                            == (start != node.identifier.coordinates))
                            || !enforce_non_match)
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
                    node.identifier.coordinates.is_secondary()
                        && ((node.identifier.has_non_match
                            == (start != node.identifier.coordinates))
                            || !enforce_non_match)
                })
                .map(|node| {
                    (
                        SecondaryAnchor::new_from_start(&node.identifier.coordinates),
                        node.cost,
                    )
                }),
        );
    }
}
