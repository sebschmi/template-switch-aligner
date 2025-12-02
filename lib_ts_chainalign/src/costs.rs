use std::ops::Range;

pub struct GapAffineCosts<Cost> {
    pub substitution: Cost,
    pub gap_open: Cost,
    pub gap_extend: Cost,
}

pub struct TsLimits {
    pub jump_12: Range<isize>,
    pub jump_34: Range<isize>,
    pub length_23: Range<usize>,
    pub ancestor_gap: Range<isize>,
}

pub struct AlignmentCosts<Cost> {
    pub primary_costs: GapAffineCosts<Cost>,
    pub secondary_costs: GapAffineCosts<Cost>,
    pub ts_base_cost: Cost,
    pub ts_limits: TsLimits,
}

impl<Cost> GapAffineCosts<Cost> {
    pub fn new(substitution: Cost, gap_open: Cost, gap_extend: Cost) -> Self {
        Self {
            substitution,
            gap_open,
            gap_extend,
        }
    }
}

impl<Cost> AlignmentCosts<Cost> {
    pub fn new(
        primary_costs: GapAffineCosts<Cost>,
        secondary_costs: GapAffineCosts<Cost>,
        ts_base_cost: Cost,
        ts_limits: TsLimits,
    ) -> Self {
        Self {
            primary_costs,
            secondary_costs,
            ts_base_cost,
            ts_limits,
        }
    }
}
