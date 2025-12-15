use compact_genome::interface::alphabet::Alphabet;
use generic_a_star::cost::U32Cost;
use lib_tsalign::a_star_aligner::{
    alignment_geometry::AlignmentRange, alignment_result::AlignmentResult,
    template_switch_distance::AlignmentType,
};
use log::{debug, info, trace};

use crate::{
    alignment::{coordinates::AlignmentCoordinates, sequences::AlignmentSequences},
    anchors::Anchors,
    chain_align::AlignmentParameters,
    chaining_cost_function::ChainingCostFunction,
    chaining_lower_bounds::ChainingLowerBounds,
    costs::AlignmentCosts,
};

pub mod alignment;
pub mod anchors;
pub mod chain_align;
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

#[expect(clippy::too_many_arguments)]
pub fn align<AlphabetType: Alphabet>(
    reference: Vec<u8>,
    query: Vec<u8>,
    range: AlignmentRange,
    parameters: &AlignmentParameters,
    rc_fn: &dyn Fn(u8) -> u8,
    reference_name: &str,
    query_name: &str,
    chaining_lower_bounds: &ChainingLowerBounds<U32Cost>,
) -> AlignmentResult<AlignmentType, U32Cost> {
    debug!(
        "Reference sequence: {}",
        String::from_utf8_lossy(&reference)
    );
    debug!("Query sequence: {}", String::from_utf8_lossy(&query));
    info!("Aligning on subsequence {}", range);

    let sequences = AlignmentSequences::new_named(
        reference,
        query,
        reference_name.to_string(),
        query_name.to_string(),
    );
    let k = chaining_lower_bounds.max_match_run() + 1;

    let anchors = Anchors::new(&sequences, range.clone(), k, rc_fn);
    trace!("Anchors:\n{anchors}");
    let start = AlignmentCoordinates::new_primary(range.reference_offset(), range.query_offset());
    let end = AlignmentCoordinates::new_primary(range.reference_limit(), range.query_limit());
    let mut chaining_cost_function =
        ChainingCostFunction::new_from_lower_bounds(chaining_lower_bounds, &anchors, start, end);

    chain_align::align::<AlphabetType, _>(
        &sequences,
        start,
        end,
        parameters,
        chaining_lower_bounds.alignment_costs(),
        rc_fn,
        chaining_lower_bounds.max_match_run(),
        &anchors,
        &mut chaining_cost_function,
    )
}
