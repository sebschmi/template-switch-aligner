use generic_a_star::cost::U32Cost;
use lib_tsalign::a_star_aligner::{
    alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    template_switch_distance::AlignmentType,
};

use crate::{chaining_lower_bounds::ChainingLowerBounds, costs::AlignmentCosts};

pub mod alignment;
pub mod chaining_cost_function;
pub mod chaining_lower_bounds;
pub mod costs;
pub mod exact_chaining;

/// Perform preprocessing for tschainalign.
///
/// * `max_n` is the maximum sequence length that the lower bounds should support.
/// * `max_match_run` is the maximum consecutive sequence of matches that is allowed.
///   Set this to `k-1`, if the anchors are `k`-mers.
/// * `alignment_costs` is the cost function for the alignment.
pub fn preprocess(
    max_n: usize,
    max_match_run: u32,
    alignment_costs: AlignmentCosts<U32Cost>,
) -> ChainingLowerBounds<U32Cost> {
    ChainingLowerBounds::new(max_n, max_match_run, alignment_costs)
}

pub fn align(
    _reference: &[u8],
    _query: &[u8],
    _range: AlignmentRange,
    _reference_name: &str,
    _query_name: &str,
    _chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<AlignmentType, U32Cost> {
    todo!()
}
