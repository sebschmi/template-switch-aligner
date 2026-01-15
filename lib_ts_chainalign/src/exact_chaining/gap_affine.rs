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

    fn allow_direct_chaining(
        &self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> bool {
        start == self.sequences.primary_start() || end == self.sequences.primary_end()
    }

    fn allow_all_matches(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        let minimum_primary_sequence_length =
            (self.sequences.primary_end().primary_ordinate_a().unwrap()
                - self.sequences.primary_start().primary_ordinate_a().unwrap())
            .min(
                self.sequences.primary_end().primary_ordinate_b().unwrap()
                    - self.sequences.primary_start().primary_ordinate_b().unwrap(),
            );
        start == self.sequences.primary_start()
            && end == self.sequences.primary_end()
            && u32::try_from(minimum_primary_sequence_length).unwrap() <= self.max_match_run
    }

    /// Align from start to end.
    ///
    /// Additionally continue the alignment to all nodes with the same cost as the alignment cost from start to end.
    ///
    /// Collect all closed nodes into the given output lists.
    /// Note that the output lists may contain duplicate anchors with different cost.
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

        let context = Context::new(
            self.cost_table,
            self.sequences,
            self.rc_fn,
            start,
            end,
            self.allow_direct_chaining(start, end),
            self.allow_all_matches(start, end),
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

        self.fill_additional_targets(
            &a_star,
            start,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());

        (cost, alignment)
    }

    /// Align from start until the cost limit is reached.
    ///
    /// Collect all closed nodes into the given output lists.
    /// Note that the output lists may contain duplicate anchors with different cost.
    pub fn align_until_cost_limit(
        &mut self,
        start: AlignmentCoordinates,
        cost_limit: Cost,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) {
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
            true,
            true,
            self.max_match_run,
        );
        let mut a_star = AStar::new_with_buffers(context, self.a_star_buffers.take().unwrap());
        a_star.initialise();
        a_star.search_until_with_target_policy(|_, node| node.cost > cost_limit, true);

        self.fill_additional_targets(
            &a_star,
            start,
            additional_primary_targets_output,
            additional_secondary_targets_output,
        );
        self.a_star_buffers = Some(a_star.into_buffers());
    }

    fn fill_additional_targets(
        &self,
        a_star: &AStar<Context<Cost>>,
        start: AlignmentCoordinates,
        additional_primary_targets_output: &mut impl Extend<(PrimaryAnchor, Cost)>,
        additional_secondary_targets_output: &mut impl Extend<(SecondaryAnchor, Cost)>,
    ) {
        additional_primary_targets_output.extend(
            a_star
                .iter_closed_nodes()
                .filter(|node| node.identifier.coordinates.is_primary())
                .filter(|node| {
                    start != node.identifier.coordinates
                        || self.allow_direct_chaining(start, node.identifier.coordinates)
                })
                .filter(|node| {
                    start == node.identifier.coordinates
                        || node.identifier.has_non_match
                        || self.allow_all_matches(start, node.identifier.coordinates)
                })
                .filter(|node| {
                    // Filter duplicates.
                    let mut pair_identifier = node.identifier;
                    pair_identifier.has_non_match = !pair_identifier.has_non_match;
                    a_star
                        .closed_node(&pair_identifier)
                        .map(|pair| {
                            pair.cost > node.cost
                                || (pair.cost == node.cost && node.identifier.has_non_match)
                        })
                        .unwrap_or(true)
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
                .filter(|node| node.identifier.coordinates.is_secondary())
                .filter(|node| {
                    start != node.identifier.coordinates
                        || self.allow_direct_chaining(start, node.identifier.coordinates)
                })
                .filter(|node| {
                    start == node.identifier.coordinates
                        || node.identifier.has_non_match
                        || self.allow_all_matches(start, node.identifier.coordinates)
                })
                .filter(|node| {
                    // Filter duplicates.
                    let mut pair_identifier = node.identifier;
                    pair_identifier.has_non_match = !pair_identifier.has_non_match;
                    a_star
                        .closed_node(&pair_identifier)
                        .map(|pair| {
                            pair.cost > node.cost
                                || (pair.cost == node.cost && node.identifier.has_non_match)
                        })
                        .unwrap_or(true)
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
