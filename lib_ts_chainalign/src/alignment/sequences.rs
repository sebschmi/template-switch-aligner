use crate::alignment::{
    coordinates::AlignmentCoordinates,
    ts_kind::{TsAncestor, TsDescendant},
};

pub struct AlignmentSequences {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
}

impl AlignmentSequences {
    pub fn new(seq1: Vec<u8>, seq2: Vec<u8>) -> Self {
        Self { seq1, seq2 }
    }

    pub fn characters(
        &self,
        coordinates: AlignmentCoordinates,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> (u8, u8) {
        match coordinates {
            AlignmentCoordinates::Primary { a, b } => (self.seq1[a], self.seq2[b]),
            AlignmentCoordinates::Secondary {
                ancestor,
                descendant,
                ts_kind,
            } => (
                match ts_kind.ancestor {
                    TsAncestor::Seq1 => self.seq1[ancestor],
                    TsAncestor::Seq2 => self.seq2[ancestor],
                },
                rc_fn(match ts_kind.descendant {
                    TsDescendant::Seq1 => self.seq1[descendant],
                    TsDescendant::Seq2 => self.seq2[descendant],
                }),
            ),
        }
    }

    pub fn start(&self) -> AlignmentCoordinates {
        AlignmentCoordinates::new_primary(0, 0)
    }

    pub fn end(&self) -> AlignmentCoordinates {
        AlignmentCoordinates::new_primary(self.seq1.len(), self.seq2.len())
    }
}
