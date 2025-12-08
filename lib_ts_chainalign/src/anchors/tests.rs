use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;

use crate::{
    alignment::sequences::AlignmentSequences,
    anchors::{Anchors, PrimaryAnchor, SecondaryAnchor},
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

#[test]
fn test_coordinates() {
    let sequences = AlignmentSequences::new(b"ACAC".to_vec(), b"ACGT".to_vec());
    let range = AlignmentRange::new_complete(sequences.seq1().len(), sequences.seq2().len());
    let k = 2;

    let anchors = Anchors::new(&sequences, range, k, &rc_fn);
    assert_eq!(anchors.primary, [(0, 0), (2, 0)].map(PrimaryAnchor::from));
    assert!(anchors.secondary_11.is_empty());
    assert_eq!(
        anchors.secondary_12,
        [(2, 2), (4, 2)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_21,
        [(4, 0), (4, 2)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_22,
        [(4, 0), (3, 1), (2, 2)].map(SecondaryAnchor::from)
    );
}

#[test]
fn test_coordinates_rev() {
    let sequences = AlignmentSequences::new(b"ACGT".to_vec(), b"ACAC".to_vec());
    let range = AlignmentRange::new_complete(sequences.seq1().len(), sequences.seq2().len());
    let k = 2;

    let anchors = Anchors::new(&sequences, range, k, &rc_fn);
    assert_eq!(anchors.primary, [(0, 0), (0, 2)].map(PrimaryAnchor::from));
    assert!(anchors.secondary_22.is_empty());
    assert_eq!(
        anchors.secondary_21,
        [(2, 2), (4, 2)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_12,
        [(4, 0), (4, 2)].map(SecondaryAnchor::from)
    );
    assert_eq!(
        anchors.secondary_11,
        [(4, 0), (3, 1), (2, 2)].map(SecondaryAnchor::from)
    );
}
