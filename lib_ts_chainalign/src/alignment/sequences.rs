use crate::alignment::coordinates::AlignmentCoordinates;

pub struct AlignmentSequences {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
}

impl AlignmentSequences {
    pub fn new(seq1: Vec<u8>, seq2: Vec<u8>) -> Self {
        Self { seq1, seq2 }
    }

    pub fn characters(&self, coordinates: AlignmentCoordinates) -> (u8, u8) {
        (
            self.seq1[coordinates.seq1().ordinate()],
            self.seq2[coordinates.seq2().ordinate()],
        )
    }
}
