use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::Cost;
use ndarray::Array2;

use crate::{
    a_star_aligner::template_switch_distance::strategies::{
        chaining::NoChainingStrategy, node_ord::CostOnlyNodeOrdStrategy,
        secondary_deletion_strategy::AllowSecondaryDeletionStrategy,
        shortcut::TemplateSwitchShortcutStrategy,
        template_switch_count::NoTemplateSwitchCountStrategy,
        template_switch_min_length::NoTemplateSwitchMinLengthStrategy, AlignmentStrategySelection,
    },
    config::TemplateSwitchConfig,
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
};

use super::template_switch::TemplateSwitchLowerBoundMatrix;

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemplateSwitchAlignmentLowerBoundMatrix {
    matrix: Array2<Cost>,
    shift_x: isize,
    shift_y: isize,
}

type TSALBAlignmentStrategies<AlphabetType> = AlignmentStrategySelection<
    AlphabetType,
    CostOnlyNodeOrdStrategy,
    NoTemplateSwitchMinLengthStrategy,
    NoChainingStrategy,
    NoTemplateSwitchCountStrategy,
    AllowSecondaryDeletionStrategy,
    TemplateSwitchShortcutStrategy,
>;

impl TemplateSwitchAlignmentLowerBoundMatrix {
    pub fn new<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
        tslb_matrix: &TemplateSwitchLowerBoundMatrix,
    ) -> Self {
        let lower_bound_config = generate_template_switch_alignment_lower_bound_config(config);
        todo!()
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
