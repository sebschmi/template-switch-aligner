use std::{collections::HashMap, iter};

use compact_genome::{
    implementation::vec_sequence::{SliceSubGenome, VectorGenome},
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{cost::Cost, AStar, AStarNode, AStarResult};
use log::{debug, info, trace};
use ndarray::Array2;

use crate::{
    a_star_aligner::{
        alignment_result::IAlignmentType,
        template_switch_distance::{
            context::Memory,
            strategies::{
                chaining::NoChainingStrategy, node_ord::CostOnlyNodeOrdStrategy,
                secondary_deletion_strategy::AllowSecondaryDeletionStrategy,
                shortcut::TemplateSwitchLowerBoundShortcutStrategy,
                template_switch_count::NoTemplateSwitchCountStrategy,
                template_switch_min_length::NoTemplateSwitchMinLengthStrategy,
                AlignmentStrategySelection,
            },
            Context, Identifier,
        },
        AlignmentContext,
    },
    config::TemplateSwitchConfig,
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
};

use super::template_switch::TemplateSwitchLowerBoundMatrix;

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchAlignmentLowerBoundMatrix {
    matrix: Array2<Cost>,
}

type TSALBAlignmentStrategies<AlphabetType> = AlignmentStrategySelection<
    AlphabetType,
    CostOnlyNodeOrdStrategy,
    NoTemplateSwitchMinLengthStrategy,
    NoChainingStrategy,
    NoTemplateSwitchCountStrategy,
    AllowSecondaryDeletionStrategy,
    TemplateSwitchLowerBoundShortcutStrategy,
>;

impl TemplateSwitchAlignmentLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
        tslb_matrix: &TemplateSwitchLowerBoundMatrix,
        reference_limit: usize,
        query_limit: usize,
    ) -> Self {
        info!("Computing TS alignment lower bound matrix...");
        let lower_bound_config = generate_template_switch_alignment_lower_bound_config(config);

        let mut closed_lower_bounds = HashMap::new();
        let genome_length = reference_limit.max(query_limit);

        debug!("Using genome length {genome_length}");
        assert!(genome_length < usize::try_from(isize::MAX).unwrap() / 2);
        let genome = VectorGenome::<AlphabetType>::from_iter(
            iter::repeat(AlphabetType::iter().next().unwrap()).take(genome_length),
        );
        let mut a_star = AStar::new(Context::<_, TSALBAlignmentStrategies<AlphabetType>>::new(
            genome.as_genome_subsequence(),
            genome.as_genome_subsequence(),
            lower_bound_config.clone(),
            Memory {
                template_switch_min_length: (),
                chaining: (),
                template_switch_count: (),
                shortcut: tslb_matrix.clone(),
            },
        ));
        a_star.initialise();

        for (target_reference_index, target_query_index) in
            (0..reference_limit).flat_map(|reference_index| {
                (0..query_limit).map(move |query_index| (reference_index, query_index))
            })
        {
            if closed_lower_bounds.contains_key(&(target_reference_index, target_query_index)) {
                continue;
            }

            trace!("Searching for target ({target_reference_index}, {target_query_index})");

            match a_star.search_until(|_, node| match *node.identifier() {
                Identifier::Primary {
                    reference_index,
                    query_index,
                    ..
                }
                | Identifier::PrimaryReentry {
                    reference_index,
                    query_index,
                    ..
                } => {
                    trace!("Found target {} at cost {}", node.identifier(), node.cost());

                    let previous =
                        closed_lower_bounds.insert((reference_index, query_index), node.cost());
                    debug_assert!(previous.is_none(), "Search finds each target at most once.");

                    reference_index == target_reference_index && query_index == target_query_index
                }
                Identifier::TemplateSwitchEntrance { .. }
                | Identifier::Secondary { .. }
                | Identifier::TemplateSwitchExit { .. } => {
                    unreachable!()
                }
            }) {
                AStarResult::FoundTarget { identifier, cost } => {
                    trace!("Search termianted with target {identifier} at cost {cost}");

                    if let Identifier::Primary { .. } | Identifier::PrimaryReentry { .. } =
                        identifier
                    {
                        debug_assert_eq!(
                            closed_lower_bounds.get(&(target_reference_index, target_query_index)),
                            Some(cost).as_ref()
                        );
                    } else {
                        trace!("{:?}", {
                            let mut alignment = Vec::new();

                            // Backtrack.
                            for alignment_type in
                                a_star.backtrack().into_iter().map(
                                    <Context<
                                        SliceSubGenome<_>,
                                        TSALBAlignmentStrategies<AlphabetType>,
                                    > as AlignmentContext>::AlignmentType::from,
                                )
                            {
                                if let Some((count, previous_alignment_type)) = alignment.last_mut()
                                {
                                    if alignment_type.is_repeated(previous_alignment_type) {
                                        *count += 1;
                                    } else {
                                        alignment.push((1, alignment_type));
                                    }
                                } else {
                                    alignment.push((1, alignment_type));
                                }
                            }

                            alignment.reverse();
                            alignment
                        });
                    }
                }
                AStarResult::NoTarget => {
                    unreachable!("Search terminated without target");
                }
            }
        }

        let (max_x, max_y) = closed_lower_bounds
            .iter()
            .fold((usize::MIN, usize::MIN), |(max_x, max_y), (&(x, y), _)| {
                (max_x.max(x), max_y.max(y))
            });

        let mut matrix = Array2::from_elem((max_x + 1, max_y + 1), Cost::MAX);
        for ((x, y), cost) in matrix.indexed_iter_mut() {
            *cost = *closed_lower_bounds.get(&(x, y)).unwrap();
            assert_ne!(*cost, Cost::MAX);
        }

        Self { matrix }
    }
}

fn generate_template_switch_alignment_lower_bound_config<AlphabetType: Alphabet>(
    config: &TemplateSwitchConfig<AlphabetType>,
) -> TemplateSwitchConfig<AlphabetType> {
    // Template switches are handled by the TS lower bound matrix, hence we disable them in the search.
    TemplateSwitchConfig {
        left_flank_length: config.left_flank_length,
        right_flank_length: config.right_flank_length,
        min_length: usize::MAX,

        base_cost: Cost::MAX,

        primary_edit_costs: config.primary_edit_costs.clone(),
        secondary_edit_costs: GapAffineAlignmentCostTable::new_max(),
        left_flank_edit_costs: config.left_flank_edit_costs.clone(),
        right_flank_edit_costs: config.right_flank_edit_costs.clone(),

        offset_costs: CostFunction::new_max(),
        length_costs: CostFunction::new_max(),
        length_difference_costs: CostFunction::new_max(),
    }
}
