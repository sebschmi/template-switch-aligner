use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U64Cost;
use num_traits::real::Real;

use super::{gap_affine_edit_distance::ScoringTable, gap_affine_edit_distance_a_star_align};

#[test]
fn match_overtakes_gap() {
    let reference = VectorGenome::<DnaAlphabet>::from_iter_u8("AGT".bytes()).unwrap();
    let query = VectorGenome::from_iter_u8("GTCC".bytes()).unwrap();
    let scoring_table = ScoringTable::<U64Cost> {
        match_cost: 0u64.into(),
        substitution_cost: 2u64.into(),
        gap_open_cost: 4u64.into(),
        gap_extend_cost: 1u64.into(),
    };

    let alignment_result = gap_affine_edit_distance_a_star_align(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        scoring_table,
    );

    assert_eq!(alignment_result.cigar(), "1D2M2I");
    assert!((alignment_result.statistics().cost - 9.0).abs() < 1e-6);
}
