use std::{
    collections::{HashMap, HashSet},
    iter,
};

use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{alphabet::Alphabet, sequence::GenomeSequence},
};
use generic_a_star::{cost::Cost, AStar};
use ndarray::Array2;

use crate::{
    a_star_aligner::template_switch_distance::{Context, GapType, Identifier, Node},
    config::TemplateSwitchConfig,
    costs::gap_affine::GapAffineAlignmentCostTable,
};

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchLowerBoundMatrix {
    matrix: Array2<Cost>,
    min_distance_between_two_template_switches: usize,
    template_switch_min_length: isize,
    template_switch_max_length: isize,
}

impl TemplateSwitchLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(config: &TemplateSwitchConfig<AlphabetType>) -> Self {
        let lower_bound_config = generate_template_switch_lower_bound_config(config);
        let mut open_lower_bounds = HashSet::new();
        open_lower_bounds.insert((0, 0));
        let mut closed_lower_bounds = HashMap::new();
        let mut string_length = 10_000;

        loop {
            let string = VectorGenome::from_iter(
                iter::repeat(AlphabetType::iter().next().unwrap()).take(string_length),
            );
            let mut a_star = AStar::new(Context::new(
                string.as_genome_subsequence(),
                string.as_genome_subsequence(),
                lower_bound_config,
            ));
            let root_x = string_length / 2;
            let root_y = string_length / 2;
            a_star.initialise_with(|context| Node::new_root_at(root_x, root_y, context));

            while let Some(coordinates) = open_lower_bounds.iter().next() {
                let (x, y) = open_lower_bounds.take(coordinates).unwrap();
                debug_assert!(!closed_lower_bounds.contains_key(&(x, y)));

                todo!("search shortest path to (x, y), while aborting in case the string is too short");

                todo!("insert neighbours into open lower bounds, if they are not yet closed");
            }
        }

        let min_distance_between_two_template_switches =
            usize::try_from(config.left_flank_length + config.right_flank_length).unwrap();

        Self {
            matrix,
            min_distance_between_two_template_switches,
            template_switch_min_length,
            template_switch_max_length,
        }
    }

    pub fn min_distance_between_two_template_switches(&self) -> usize {
        self.min_distance_between_two_template_switches
    }

    pub fn template_switch_min_length(&self) -> isize {
        self.template_switch_min_length
    }

    pub fn template_switch_max_length(&self) -> isize {
        self.template_switch_max_length
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
