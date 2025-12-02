use generic_a_star::cost::U32Cost;

use crate::{
    alignment::{
        AlignmentType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
        ts_kind::TsKind,
    },
    costs::{AlignmentCosts, GapAffineCosts, TsLimits},
    exact_chaining::ts_12_jump::Ts12JumpAlignment,
};

fn rc_fn(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        c => unimplemented!("Unsupported character {c}"),
    }
}

#[ignore]
#[test]
fn test_start_end() {
    let seq1 = b"AAGT".to_vec();
    let seq2 = b"ACGTT".to_vec();
    let sequences = AlignmentSequences::new(seq1, seq2);
    let cost_table = AlignmentCosts {
        primary_costs: GapAffineCosts::new(
            U32Cost::from(2u8),
            U32Cost::from(3u8),
            U32Cost::from(1u8),
        ),
        secondary_costs: GapAffineCosts::new(
            U32Cost::from(4u8),
            U32Cost::from(6u8),
            U32Cost::from(2u8),
        ),
        ts_base_cost: U32Cost::from(2u8),
        ts_limits: TsLimits {
            jump_12: -100..100,
            jump_34: -100..100,
            length_23: 0..100,
            ancestor_gap: -100..100,
        },
    };

    let start = AlignmentCoordinates::new_primary(0, 0);
    let end = AlignmentCoordinates::new_secondary(4, 0, TsKind::TS21);
    let alignment = Ts12JumpAlignment::new(start, end, &sequences, &cost_table, &rc_fn, u32::MAX);

    assert_eq!(alignment.start(), start);
    assert_eq!(alignment.end(), end);
    assert_eq!(
        alignment.alignment().alignment,
        vec![
            (1, AlignmentType::Match),
            (1, AlignmentType::Substitution),
            (1, AlignmentType::Match),
            (
                1,
                AlignmentType::TsStart {
                    jump: -1,
                    ts_kind: TsKind::TS21
                }
            ),
            (2, AlignmentType::Match),
        ]
    );
    assert_eq!(alignment.cost(), U32Cost::from(4u8));
}
