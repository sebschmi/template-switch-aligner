use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;

use super::MultipairAlignmentRenderer;

#[test]
fn test_parallel_gaps() {
    let mut renderer = MultipairAlignmentRenderer::new("B".to_string(), "GGG");
    renderer.add_aligned_sequence(
        &"B".to_string(),
        0,
        "A".to_string(),
        "GGGG",
        [
            (2, AlignmentType::PrimaryMatch),
            (1, AlignmentType::PrimaryInsertion),
            (1, AlignmentType::PrimaryMatch),
        ],
        true,
        false,
    );
    renderer.add_aligned_sequence(
        &"B".to_string(),
        0,
        "C".to_string(),
        "GGGG",
        [
            (2, AlignmentType::PrimaryMatch),
            (1, AlignmentType::PrimaryInsertion),
            (1, AlignmentType::PrimaryMatch),
        ],
        true,
        false,
    );

    let mut output = Vec::new();
    renderer
        .render(
            &mut output,
            [&"A".to_string(), &"B".to_string(), &"C".to_string()],
        )
        .unwrap();

    let output = String::from_utf8(output).unwrap();
    assert_eq!(output.trim(), "A: GGGG\nB: GG-G\nC: GGGG");
}
