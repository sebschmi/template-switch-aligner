use crate::alignment::{
    coordinates::AlignmentCoordinates,
    ts_kind::{TsAncestor, TsDescendant, TsKind},
};

pub struct AlignmentSequences {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
    seq1_name: String,
    seq2_name: String,
}

impl AlignmentSequences {
    pub fn new(seq1: Vec<u8>, seq2: Vec<u8>) -> Self {
        Self {
            seq1,
            seq2,
            seq1_name: "seq1".to_string(),
            seq2_name: "seq2".to_string(),
        }
    }

    pub fn new_named(seq1: Vec<u8>, seq2: Vec<u8>, seq1_name: String, seq2_name: String) -> Self {
        Self {
            seq1,
            seq2,
            seq1_name,
            seq2_name,
        }
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
                    TsAncestor::Seq1 => self.seq1[ancestor - 1],
                    TsAncestor::Seq2 => self.seq2[ancestor - 1],
                },
                rc_fn(match ts_kind.descendant {
                    TsDescendant::Seq1 => self.seq1[descendant],
                    TsDescendant::Seq2 => self.seq2[descendant],
                }),
            ),
        }
    }

    pub fn primary_start(&self) -> AlignmentCoordinates {
        AlignmentCoordinates::new_primary(0, 0)
    }

    pub fn primary_end(&self) -> AlignmentCoordinates {
        self.end(None)
    }

    pub fn secondary_end(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        self.end(Some(ts_kind))
    }

    pub fn end(&self, ts_kind: Option<TsKind>) -> AlignmentCoordinates {
        match ts_kind {
            None => AlignmentCoordinates::new_primary(self.seq1.len(), self.seq2.len()),
            Some(ts_kind @ (TsKind::TS11 | TsKind::TS21)) => {
                AlignmentCoordinates::new_secondary(0, self.seq1.len(), ts_kind)
            }
            Some(ts_kind @ (TsKind::TS12 | TsKind::TS22)) => {
                AlignmentCoordinates::new_secondary(0, self.seq2.len(), ts_kind)
            }
        }
    }

    pub fn seq1(&self) -> &[u8] {
        &self.seq1
    }

    pub fn seq2(&self) -> &[u8] {
        &self.seq2
    }

    pub fn seq1_name(&self) -> &str {
        &self.seq1_name
    }

    pub fn seq2_name(&self) -> &str {
        &self.seq2_name
    }
}
