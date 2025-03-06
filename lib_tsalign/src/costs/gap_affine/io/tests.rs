use compact_genome::implementation::alphabets::dna_alphabet::DnaAlphabet;
use generic_a_star::cost::U64Cost;

use crate::costs::gap_affine::GapAffineAlignmentCostTable;

#[test]
fn simple_example() {
    let input = "# Simple Example\n\nSubstitutionCostTable\n  |  A  C  G  T\n--+------------\nA | 10  9 11  7\nC |  2  1  3  4\nG |  5  6 17  8\nT | 99  0 50 51\n\nGapOpenCostVector\n A C G T\n 3 4 5 1\n\nGapExtendCostVector\n  A  C  G  T\n 10 15  0  1\n";
    let expected_parsing_result = GapAffineAlignmentCostTable::<DnaAlphabet, U64Cost> {
        name: "Simple Example".to_string(),
        substitution_cost_table: [10u64, 9, 11, 7, 2, 1, 3, 4, 5, 6, 17, 8, 99, 0, 50, 51]
            .into_iter()
            .map(Into::into)
            .collect(),
        gap_open_cost_vector: [3u64, 4, 5, 1].into_iter().map(Into::into).collect(),
        gap_extend_cost_vector: [10u64, 15, 0, 1].into_iter().map(Into::into).collect(),
        phantom_data: Default::default(),
    };

    let actual_parsing_result =
        GapAffineAlignmentCostTable::<DnaAlphabet, U64Cost>::read_plain(input.as_bytes()).unwrap();
    let mut writer = Vec::new();
    actual_parsing_result.write_plain(&mut writer).unwrap();
    let output = String::from_utf8(writer).unwrap();

    assert_eq!(expected_parsing_result, actual_parsing_result);
    assert_eq!(input, output);
}
