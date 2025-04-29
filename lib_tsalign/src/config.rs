use compact_genome::interface::alphabet::Alphabet;
use num_traits::{Bounded, bounds::UpperBounded};

use crate::{
    costs::{cost_function::CostFunction, gap_affine::GapAffineAlignmentCostTable},
    error::{Error, Result},
};

pub mod io;

#[derive(Debug, Eq, PartialEq)]
pub struct TemplateSwitchConfig<AlphabetType, Cost> {
    // Limits
    pub left_flank_length: isize,
    pub right_flank_length: isize,
    pub min_length: usize,

    // Base cost
    pub base_cost: BaseCost<Cost>,

    // Edit costs
    pub primary_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub left_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub right_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,

    // Jump costs
    pub offset_costs: CostFunction<isize, Cost>,
    pub length_costs: CostFunction<usize, Cost>,
    pub length_difference_costs: CostFunction<isize, Cost>,
    pub anti_primary_gap_costs: CostFunction<isize, Cost>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCost<Cost> {
    pub rr: Cost,
    pub rq: Cost,
    pub qr: Cost,
    pub qq: Cost,
}

impl<AlphabetType, Cost: Bounded + Ord> TemplateSwitchConfig<AlphabetType, Cost> {
    /// Returns an error if any cost function is malformed.
    pub fn verify(&self) -> Result<()> {
        if !self.offset_costs.is_v_shaped() {
            Err(Error::OffsetCostsNotVShaped)
        } else if !self.length_difference_costs.is_v_shaped() {
            Err(Error::LengthDifferenceCostsNotVShaped)
        } else {
            Ok(())
        }
    }
}

impl<Cost: UpperBounded> BaseCost<Cost> {
    pub fn new_max() -> Self {
        Self {
            rr: Cost::max_value(),
            rq: Cost::max_value(),
            qr: Cost::max_value(),
            qq: Cost::max_value(),
        }
    }
}

impl<AlphabetType: Alphabet, Cost: Clone> Clone for TemplateSwitchConfig<AlphabetType, Cost> {
    fn clone(&self) -> Self {
        Self {
            left_flank_length: self.left_flank_length,
            right_flank_length: self.right_flank_length,
            min_length: self.min_length,
            base_cost: self.base_cost.clone(),
            primary_edit_costs: self.primary_edit_costs.clone(),
            secondary_edit_costs: self.secondary_edit_costs.clone(),
            left_flank_edit_costs: self.left_flank_edit_costs.clone(),
            right_flank_edit_costs: self.right_flank_edit_costs.clone(),
            offset_costs: self.offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
            anti_primary_gap_costs: self.anti_primary_gap_costs.clone(),
        }
    }
}
