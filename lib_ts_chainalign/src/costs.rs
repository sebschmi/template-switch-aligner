use std::ops::Range;

pub struct GapAffineCosts<Cost> {
    pub substitution: Cost,
    pub gap_open: Cost,
    pub gap_extend: Cost,
}

#[expect(dead_code)]
pub struct TsLimits {
    pub jump_12: Range<isize>,
    pub jump_34: Range<isize>,
    pub length_23: Range<usize>,
}

pub struct AlignmentCosts<Cost> {
    pub primary_costs: GapAffineCosts<Cost>,
    pub secondary_costs: GapAffineCosts<Cost>,
    pub ts_base_cost: Cost,
    #[expect(dead_code)]
    pub ts_limits: TsLimits,
}
