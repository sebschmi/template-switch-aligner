use generic_a_star::cost::U32Cost;
use lib_tsalign::a_star_aligner::{
    alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    template_switch_distance::AlignmentType,
};
use log::{debug, info};

use crate::{
    alignment::sequences::AlignmentSequences, anchors::Anchors,
    chaining_lower_bounds::ChainingLowerBounds, costs::AlignmentCosts,
};

pub mod alignment;
pub mod anchors;
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
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    rc_fn: &dyn Fn(u8) -> u8,
    _reference_name: &str,
    _query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<AlignmentType, U32Cost> {
    debug!(
        "Reference sequence: {}",
        String::from_utf8_lossy(&reference)
    );
    debug!("Query sequence: {}", String::from_utf8_lossy(&query));
    info!("Aligning on subsequence {}", range);

    let sequences = AlignmentSequences::new(reference, query);
    let k = chaining_lower_bounds.max_match_run() + 1;

    let anchors = Anchors::new(&sequences, range, k, rc_fn);
    println!("Anchors:\n{anchors}");
    todo!()
}
