use generic_a_star::cost::U32Cost;

use crate::alignment::AlignmentType;
use crate::gap_filling::gap_affine::{AlignmentCoordinates, GapAffineAlignment};
use crate::{alignment::sequences::AlignmentSequences, costs::GapAffineCosts};

#[test]
fn test_start_end() {
    let seq1 = b"ACGT".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(0, 0);
    let end = AlignmentCoordinates::new_forwards(4, 5);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![(4, AlignmentType::Match), (1, AlignmentType::Gap1)]
    );
    assert_eq!(alignment.cost(), U32Cost::from(3u8));
}

#[test]
fn test_partial_alignment() {
    let seq1 = b"ACCGT".to_vec();
    let seq2 = b"ACGGTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(1, 1);
    let end = AlignmentCoordinates::new_forwards(4, 4);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match)
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(2u8));
}

#[test]
fn test_gap_directions() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(1, 1);
    let end = AlignmentCoordinates::new_forwards(11, 11);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (2, AlignmentType::Match),
            (2, AlignmentType::Gap2),
            (5, AlignmentType::Match),
            (2, AlignmentType::Gap1),
            (1, AlignmentType::Match)
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(8u8));
}

#[test]
fn test_extremity_gaps() {
    let seq1 = b"ACGCCGTGTTCT".to_vec();
    let seq2 = b"ACGGTGTTAACT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(3, 3);
    let end = AlignmentCoordinates::new_forwards(10, 10);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (2, AlignmentType::Gap2),
            (5, AlignmentType::Match),
            (2, AlignmentType::Gap1),
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(8u8));
}

#[test]
fn test_extremity_substitutions() {
    let seq1 = b"AGGGA".to_vec();
    let seq2 = b"TGGGT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(0, 0);
    let end = AlignmentCoordinates::new_forwards(5, 5);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (1, AlignmentType::Substitution),
            (3, AlignmentType::Match),
            (1, AlignmentType::Substitution),
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(4u8));
}

#[test]
fn test_substitutions_as_gaps() {
    let seq1 = b"AAAAAAAAAAAAAAAAAAAA".to_vec();
    let seq2 = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(3u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_forwards(0, 0);
    let end = AlignmentCoordinates::new_forwards(20, 20);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert!(
        alignment.alignment().alignment
            == vec![(20, AlignmentType::Gap1), (20, AlignmentType::Gap2),]
            || alignment.alignment().alignment
                == vec![(20, AlignmentType::Gap2), (20, AlignmentType::Gap1),]
    );
    assert_eq!(alignment.cost(), U32Cost::from(44u8));
}
