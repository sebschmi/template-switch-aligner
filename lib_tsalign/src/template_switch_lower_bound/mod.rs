use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    iter,
};

use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{cost::Cost, AStar, AStarNode, AStarResult};
use ndarray::Array2;

use crate::{
    a_star_aligner::template_switch_distance::{
        strategies::SimpleAlignmentStrategies, Context, Identifier, Node,
    },
    config::TemplateSwitchConfig,
    costs::gap_affine::GapAffineAlignmentCostTable,
};

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchLowerBoundMatrix {
    matrix: Array2<Cost>,
    min_distance_between_two_template_switches: usize,
}

impl TemplateSwitchLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(config: &TemplateSwitchConfig<AlphabetType>) -> Self {
        let lower_bound_config = generate_template_switch_lower_bound_config(config);
        let mut open_lower_bounds = HashSet::new();
        open_lower_bounds.insert((0isize, 0isize, false));
        let mut closed_lower_bounds = HashMap::new();
        let mut genome_length = 10_000;

        'outer: loop {
            assert!(genome_length < usize::try_from(isize::MAX).unwrap() / 2);
            let genome = VectorGenome::<AlphabetType>::from_iter(
                iter::repeat(AlphabetType::iter().next().unwrap()).take(genome_length),
            );
            let mut a_star =
                AStar::new(Context::<_, SimpleAlignmentStrategies<AlphabetType>>::new(
                    genome.as_genome_subsequence(),
                    genome.as_genome_subsequence(),
                    lower_bound_config.clone(),
                ));
            let root_xy = genome_length / 2;
            a_star.initialise_with(|context| Node::new_root_at(root_xy, root_xy, context));

            let root_xy_isize = isize::try_from(root_xy).unwrap();
            let string_length_isize = isize::try_from(genome_length).unwrap();

            while let Some(coordinates) = open_lower_bounds.iter().next().copied() {
                let (x, y, from_target) = open_lower_bounds.take(&coordinates).unwrap();
                debug_assert!(!closed_lower_bounds.contains_key(&(x, y)));
                let reference_target = usize::try_from(x + root_xy_isize).unwrap();
                let query_target = usize::try_from(y + root_xy_isize).unwrap();

                let has_target =
                    match a_star.search_until(|context, node| match *node.identifier() {
                        Identifier::Primary {
                            reference_index,
                            query_index,
                            ..
                        } => {
                            reference_index == 0
                                || query_index == 0
                                || reference_index == genome_length - 1
                                || query_index == genome_length - 1
                        }
                        Identifier::PrimaryReentry {
                            reference_index,
                            query_index,
                            ..
                        } => reference_index == reference_target && query_index == query_target,
                        Identifier::TemplateSwitchEntrance {
                            entrance_reference_index,
                            entrance_query_index,
                            template_switch_first_offset,
                            ..
                        } => {
                            entrance_reference_index == 0
                                || entrance_query_index == 0
                                || entrance_reference_index == genome_length - 1
                                || entrance_query_index == genome_length - 1
                                || template_switch_first_offset + root_xy_isize == 0
                                || template_switch_first_offset + root_xy_isize
                                    == string_length_isize - 1
                        }
                        Identifier::Secondary {
                            primary_index,
                            secondary_index,
                            ..
                        } => {
                            primary_index == 0
                                || primary_index == genome_length - 1
                                || secondary_index == 0
                                || secondary_index == genome_length - 1
                        }
                        Identifier::TemplateSwitchExit { .. } => {
                            node.generate_primary_reentry_successor(context).is_none()
                        }
                    }) {
                        AStarResult::FoundTarget { identifier, cost } => {
                            if let Identifier::PrimaryReentry { .. } = identifier {
                                let previous = closed_lower_bounds.insert((x, y), cost);
                                debug_assert!(previous.is_none());
                                true
                            } else {
                                open_lower_bounds.insert((x, y, from_target));
                                genome_length *= 2;
                                continue 'outer;
                            }
                        }
                        AStarResult::NoTarget => {
                            let previous = closed_lower_bounds.insert((x, y), Cost::MAX);
                            debug_assert!(previous.is_none());
                            false
                        }
                    };

                if has_target || !from_target {
                    let from_target = has_target;

                    for (x, y) in [
                        (x + 1, y + 1),
                        (x + 1, y),
                        (x + 1, y - 1),
                        (x, y + 1),
                        (x, y - 1),
                        (x - 1, y + 1),
                        (x - 1, y),
                        (x - 1, y - 1),
                    ] {
                        if !closed_lower_bounds.contains_key(&(x, y)) {
                            open_lower_bounds.insert((x, y, from_target));
                        }
                    }
                }
            }

            break;
        }

        let min_distance_between_two_template_switches =
            usize::try_from(config.left_flank_length + config.right_flank_length).unwrap();

        let (min_x, max_x, min_y, max_y) = closed_lower_bounds.keys().fold(
            (isize::MAX, isize::MIN, isize::MAX, isize::MIN),
            |(min_x, max_x, min_y, max_y), &(x, y)| {
                (min_x.min(x), max_x.max(x), min_y.min(y), max_y.max(y))
            },
        );

        let shift_x = -min_x;
        let shift_y = -max_y;
        let len_x = usize::try_from(max_x - min_x + 1).unwrap();
        let len_y = usize::try_from(max_y - min_y + 1).unwrap();

        let mut matrix = Array2::from_elem([len_x, len_y], Cost::MAX);
        for ((x, y), cost) in closed_lower_bounds {
            let x = usize::try_from(x + shift_x).unwrap();
            let y = usize::try_from(y + shift_y).unwrap();
            matrix[(x, y)] = cost;
        }

        Self {
            matrix,
            min_distance_between_two_template_switches,
        }
    }

    pub fn min_distance_between_two_template_switches(&self) -> usize {
        self.min_distance_between_two_template_switches
    }
}

fn generate_template_switch_lower_bound_config<AlphabetType: Alphabet>(
    config: &TemplateSwitchConfig<AlphabetType>,
) -> TemplateSwitchConfig<AlphabetType> {
    TemplateSwitchConfig {
        left_flank_length: 0,
        right_flank_length: 0,
        min_length: config.min_length,

        base_cost: config.base_cost,

        primary_edit_costs: GapAffineAlignmentCostTable::new_max(),
        secondary_edit_costs: config.secondary_edit_costs.clone().into_lower_bound(),
        left_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),
        right_flank_edit_costs: GapAffineAlignmentCostTable::new_max(),

        offset_costs: config.offset_costs.clone(),
        length_costs: config.length_costs.clone(),
        length_difference_costs: config.length_difference_costs.clone(),
    }
}

impl Display for TemplateSwitchLowerBoundMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "TemplateSwitchLowerBoundMatrix:")?;
        writeln!(
            f,
            "min distance: {}",
            self.min_distance_between_two_template_switches
        )?;
        writeln!(f, "{}", self.matrix)?;

        Ok(())
    }
}
