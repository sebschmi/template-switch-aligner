use compact_genome::{
    implementation::vec_sequence::VectorGenome,
    interface::{
        alphabet::dna_alphabet::DnaAlphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
};

use super::{gap_affine_edit_distance::ScoringTable, gap_affine_edit_distance_a_star_align};

#[test]
fn match_overtakes_gap() {
    let reference = VectorGenome::<DnaAlphabet>::from_iter_u8("AGT".bytes()).unwrap();
    let query = VectorGenome::from_iter_u8("GTCC".bytes()).unwrap();
    let scoring_table = ScoringTable {
        match_cost: 0.into(),
        substitution_cost: 2.into(),
        gap_open_cost: 4.into(),
        gap_extend_cost: 1.into(),
    };

    let alignment_result = gap_affine_edit_distance_a_star_align(
        reference.as_genome_subsequence(),
        query.as_genome_subsequence(),
        scoring_table,
    );

    assert_eq!(alignment_result.cigar(), "1D2M2I");
    assert_eq!(alignment_result.cost, 9.into());
}