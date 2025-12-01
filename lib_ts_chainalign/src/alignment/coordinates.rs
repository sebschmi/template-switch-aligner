use std::{fmt::Display, ops::Range};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct AlignmentCoordinates {
    seq1: SequenceOrdinate,
    seq2: SequenceOrdinate,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct SequenceOrdinate {
    ordinate: usize,
    rc: bool,
}

impl AlignmentCoordinates {
    pub fn new_forwards(seq1: usize, seq2: usize) -> Self {
        Self {
            seq1: SequenceOrdinate::new(seq1, false),
            seq2: SequenceOrdinate::new(seq2, false),
        }
    }

    pub fn seq1(&self) -> SequenceOrdinate {
        self.seq1
    }

    pub fn seq2(&self) -> SequenceOrdinate {
        self.seq2
    }

    pub fn can_increment_1(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        self.seq1()
            .can_increment(start.seq1().ordinate()..end.seq1().ordinate())
    }

    pub fn can_increment_2(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        self.seq2()
            .can_increment(start.seq2().ordinate()..end.seq2().ordinate())
    }

    pub fn can_increment_both(
        &self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> bool {
        self.can_increment_1(start, end) && self.can_increment_2(start, end)
    }

    pub fn increment_1(&self) -> Self {
        Self {
            seq1: SequenceOrdinate::new(
                if self.seq1.is_rc() {
                    self.seq1.ordinate().wrapping_sub(1)
                } else {
                    self.seq1.ordinate() + 1
                },
                self.seq1.is_rc(),
            ),
            seq2: self.seq2,
        }
    }

    pub fn increment_2(&self) -> Self {
        Self {
            seq1: self.seq1,
            seq2: SequenceOrdinate::new(
                if self.seq2.is_rc() {
                    self.seq2.ordinate().wrapping_sub(1)
                } else {
                    self.seq2.ordinate() + 1
                },
                self.seq2.is_rc(),
            ),
        }
    }

    pub fn increment_both(&self) -> Self {
        self.increment_1().increment_2()
    }
}

impl SequenceOrdinate {
    pub fn new(ordinate: usize, rc: bool) -> Self {
        Self { ordinate, rc }
    }

    /// Returns the sequence index of this ordinate.
    ///
    /// If the index runs over the end of the sequence, it may roll around to `usize::MAX` in reverse complements.
    pub fn ordinate(&self) -> usize {
        self.ordinate
    }

    pub fn is_rc(&self) -> bool {
        self.rc
    }

    pub fn can_increment(&self, range: Range<usize>) -> bool {
        range.contains(&self.ordinate)
    }
}

impl Display for SequenceOrdinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.ordinate, if self.rc { "rc" } else { "" })
    }
}
