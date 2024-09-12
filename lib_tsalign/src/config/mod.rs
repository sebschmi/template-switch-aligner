use crate::costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable};

pub mod io;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct TemplateSwitchConfig<AlphabetType> {
    // Limits
    pub left_flank_length: isize,
    pub right_flank_length: isize,

    // Edit costs
    pub primary_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub secondary_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub left_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,
    pub right_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType>,

    // Jump costs
    pub offset1_costs: CostFunction<isize>,
}
