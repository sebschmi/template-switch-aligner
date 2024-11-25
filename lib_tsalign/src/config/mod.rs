use compact_genome::interface::alphabet::Alphabet;

use crate::costs::{
    cost::Cost, cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable,
};

pub mod io;

#[derive(Debug, Eq, PartialEq)]
pub struct TemplateSwitchConfig<AlphabetType> {
    // Limits
    pub left_flank_length: isize,
    pub right_flank_length: isize,
    pub min_length: usize,

    // Base cost
    pub base_cost: Cost,

    // Edit costs
    pub primary_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub secondary_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub left_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub right_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,

    // Jump costs
    pub offset_costs: CostFunction<isize>,
    pub length_costs: CostFunction<usize>,
    pub length_difference_costs: CostFunction<isize>,
}

impl<AlphabetType: Alphabet> Clone for TemplateSwitchConfig<AlphabetType> {
    fn clone(&self) -> Self {
        Self {
            left_flank_length: self.left_flank_length,
            right_flank_length: self.right_flank_length,
            min_length: self.min_length,
            base_cost: self.base_cost,
            primary_edit_costs: self.primary_edit_costs.clone(),
            secondary_edit_costs: self.secondary_edit_costs.clone(),
            left_flank_edit_costs: self.left_flank_edit_costs.clone(),
            right_flank_edit_costs: self.right_flank_edit_costs.clone(),
            offset_costs: self.offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
        }
    }
}
