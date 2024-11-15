use crate::costs::{
    cost::Cost, cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable,
};

pub mod io;

#[derive(Debug, Clone, Eq, PartialEq)]
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
