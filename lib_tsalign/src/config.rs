use compact_genome::interface::alphabet::Alphabet;
use num_traits::{Bounded, bounds::UpperBounded};

use crate::{
    a_star_aligner::template_switch_distance::TemplateSwitchDirection,
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
    pub secondary_forward_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub secondary_reverse_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub left_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,
    pub right_flank_edit_costs: GapAffineAlignmentCostTable<AlphabetType, Cost>,

    // Jump costs
    pub offset_costs: CostFunction<isize, Cost>,
    pub length_costs: CostFunction<usize, Cost>,
    pub length_difference_costs: CostFunction<isize, Cost>,
    pub forward_anti_primary_gap_costs: CostFunction<isize, Cost>,
    pub reverse_anti_primary_gap_costs: CostFunction<isize, Cost>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCost<Cost> {
    /// Primary: reference; secondary: reference; direction: forward.
    pub rrf: Cost,
    /// Primary: reference; secondary: query; direction: forward.
    pub rqf: Cost,
    /// Primary: query; secondary: reference; direction: forward.
    pub qrf: Cost,
    /// Primary: query; secondary: query; direction: forward.
    pub qqf: Cost,
    /// Primary: reference; secondary: reference; direction: reverse.
    pub rrr: Cost,
    /// Primary: reference; secondary: query; direction: reverse.
    pub rqr: Cost,
    /// Primary: query; secondary: reference; direction: reverse.
    pub qrr: Cost,
    /// Primary: query; secondary: query; direction: reverse.
    pub qqr: Cost,
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

impl<AlphabetType, Cost> TemplateSwitchConfig<AlphabetType, Cost> {
    pub fn secondary_edit_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &GapAffineAlignmentCostTable<AlphabetType, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.secondary_forward_edit_costs,
            TemplateSwitchDirection::Reverse => &self.secondary_reverse_edit_costs,
        }
    }

    pub fn anti_primary_gap_costs(
        &self,
        direction: TemplateSwitchDirection,
    ) -> &CostFunction<isize, Cost> {
        match direction {
            TemplateSwitchDirection::Forward => &self.forward_anti_primary_gap_costs,
            TemplateSwitchDirection::Reverse => &self.reverse_anti_primary_gap_costs,
        }
    }
}

impl<Cost: UpperBounded> BaseCost<Cost> {
    pub fn new_max() -> Self {
        Self {
            rrf: Cost::max_value(),
            rqf: Cost::max_value(),
            qrf: Cost::max_value(),
            qqf: Cost::max_value(),
            rrr: Cost::max_value(),
            rqr: Cost::max_value(),
            qrr: Cost::max_value(),
            qqr: Cost::max_value(),
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
            secondary_forward_edit_costs: self.secondary_forward_edit_costs.clone(),
            secondary_reverse_edit_costs: self.secondary_reverse_edit_costs.clone(),
            left_flank_edit_costs: self.left_flank_edit_costs.clone(),
            right_flank_edit_costs: self.right_flank_edit_costs.clone(),
            offset_costs: self.offset_costs.clone(),
            length_costs: self.length_costs.clone(),
            length_difference_costs: self.length_difference_costs.clone(),
            forward_anti_primary_gap_costs: self.forward_anti_primary_gap_costs.clone(),
            reverse_anti_primary_gap_costs: self.reverse_anti_primary_gap_costs.clone(),
        }
    }
}
