use generic_a_star::cost::U32Cost;

use crate::alignment::AlignmentType;
use crate::alignment::ts_kind::TsKind;
use crate::exact_chaining::gap_affine::{AlignmentCoordinates, GapAffineAlignment};
use crate::{alignment::sequences::AlignmentSequences, costs::GapAffineCosts};

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => unimplemented!("Unsupported character {c}"),
    }
}

#[test]
fn test_start_end() {
    let seq1 = b"ACGT".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(4, 5);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![(4, AlignmentType::Match), (1, AlignmentType::GapA)]
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

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(4, 4);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

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

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(11, 11);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (2, AlignmentType::Match),
            (2, AlignmentType::GapB),
            (5, AlignmentType::Match),
            (2, AlignmentType::GapA),
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

    let start = AlignmentCoordinates::new_primary(3, 3);
    let end = AlignmentCoordinates::new_primary(10, 10);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (2, AlignmentType::GapB),
            (5, AlignmentType::Match),
            (2, AlignmentType::GapA),
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

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(5, 5);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

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

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_primary(20, 20);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert!(
        alignment.alignment().alignment
            == vec![(20, AlignmentType::GapA), (20, AlignmentType::GapB),]
            || alignment.alignment().alignment
                == vec![(20, AlignmentType::GapB), (20, AlignmentType::GapA),]
    );
    assert_eq!(alignment.cost(), U32Cost::from(44u8));
}

#[test]
fn test_max_match_run_0() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, 0);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert!(
        alignment.alignment().alignment
            == vec![(8, AlignmentType::GapA), (8, AlignmentType::GapB),]
            || alignment.alignment().alignment
                == vec![(8, AlignmentType::GapB), (8, AlignmentType::GapA),]
    );
    assert_eq!(alignment.cost(), U32Cost::from(20u8));
}

#[test]
fn test_max_match_run_1() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, 1);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::GapB),
            (1, AlignmentType::Match),
            (2, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (1, AlignmentType::GapA),
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(12u8));
}

#[test]
fn test_max_match_run_2() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"AACAACCAAA".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_primary(1, 1);
    let end = AlignmentCoordinates::new_primary(9, 9);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, 2);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (2, AlignmentType::Match),
            (2, AlignmentType::Substitution),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(6u8));
}

#[test]
fn test_secondary_12() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"TTTTTTTTTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS12);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS12);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![(8, AlignmentType::Match),]
    );
    assert_eq!(alignment.cost(), U32Cost::from(0u8));
}

#[test]
fn test_secondary_21() {
    let seq1 = b"AAAAAAAAAA".to_vec();
    let seq2 = b"TTTTTTTTTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table =
        GapAffineCosts::new(U32Cost::from(2u8), U32Cost::from(3u8), U32Cost::from(1u8));

    let start = AlignmentCoordinates::new_secondary(9, 1, TsKind::TS21);
    let end = AlignmentCoordinates::new_secondary(1, 9, TsKind::TS21);
    let alignment = GapAffineAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![(8, AlignmentType::Match),]
    );
    assert_eq!(alignment.cost(), U32Cost::from(0u8));
}
