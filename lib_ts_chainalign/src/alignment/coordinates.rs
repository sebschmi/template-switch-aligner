use std::fmt::Display;

use crate::alignment::ts_kind::{TsAncestor, TsDescendant, TsKind};

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum AlignmentCoordinates {
    Primary {
        a: usize,
        b: usize,
    },
    Secondary {
        ancestor: usize,
        descendant: usize,
        ts_kind: TsKind,
    },
}

impl AlignmentCoordinates {
    pub fn new_primary(a: usize, b: usize) -> Self {
        Self::Primary { a, b }
    }

    pub fn new_secondary(ancestor: usize, descendant: usize, ts_kind: TsKind) -> Self {
        Self::Secondary {
            ancestor,
            descendant,
            ts_kind,
        }
    }

    pub fn primary_ordinate_a(&self) -> Option<usize> {
        match self {
            AlignmentCoordinates::Primary { a, .. } => Some(*a),
            AlignmentCoordinates::Secondary { .. } => None,
        }
    }

    pub fn primary_ordinate_b(&self) -> Option<usize> {
        match self {
            AlignmentCoordinates::Primary { b, .. } => Some(*b),
            AlignmentCoordinates::Secondary { .. } => None,
        }
    }

    pub fn secondary_ordinate_ancestor(&self) -> Option<usize> {
        match self {
            AlignmentCoordinates::Secondary { ancestor, .. } => Some(*ancestor),
            AlignmentCoordinates::Primary { .. } => None,
        }
    }

    pub fn secondary_ordinate_descendant(&self) -> Option<usize> {
        match self {
            AlignmentCoordinates::Secondary { descendant, .. } => Some(*descendant),
            AlignmentCoordinates::Primary { .. } => None,
        }
    }

    pub fn ts_kind(&self) -> Option<TsKind> {
        match self {
            AlignmentCoordinates::Secondary { ts_kind, .. } => Some(*ts_kind),
            AlignmentCoordinates::Primary { .. } => None,
        }
    }

    pub fn is_primary(&self) -> bool {
        matches!(self, AlignmentCoordinates::Primary { .. })
    }

    pub fn is_secondary(&self) -> bool {
        matches!(self, AlignmentCoordinates::Secondary { .. })
    }

    /// Checks if ordinate a can be incremented.
    /// In secondary alignments, ordinate a is the ancestor.
    ///
    /// If ordinate a and the boundary have the same TS kind (possibly `None`), then the check is performed normally.
    /// If the TS kinds differ, then there is a jump before the boundary, and ordinate a can always be incremented.
    pub fn can_increment_a(&self, end: AlignmentCoordinates) -> bool {
        if self.ts_kind() == end.ts_kind() {
            match self {
                AlignmentCoordinates::Primary { a, .. } => {
                    // Incrementing primary ordinate a is always a plus operation, so we only need to check the upper bound.
                    (..end.primary_ordinate_a().unwrap()).contains(a)
                }
                AlignmentCoordinates::Secondary { ancestor, .. } => {
                    // Incrementing the secondary ancestor is always a minus operation, so we only need to check the lower bound.
                    (end.secondary_ordinate_ancestor().unwrap()..).contains(ancestor)
                }
            }
        } else {
            true
        }
    }

    /// Checks if ordinate b can be incremented.
    /// In secondary alignments, ordinate b is the descendant.
    ///
    /// If ordinate b and the `end` boundary have the same TS kind (possibly `None`), then the check is performed normally.
    /// If the TS kinds differ, then there is a jump before the `end` boundary, and ordinate b can always be incremented.
    pub fn can_increment_b(&self, end: AlignmentCoordinates) -> bool {
        if self.ts_kind() == end.ts_kind() {
            // Incrementing ordinate b is always a plus operation, so we only need to check the upper bound.
            match self {
                AlignmentCoordinates::Primary { b, .. } => {
                    (..end.primary_ordinate_b().unwrap()).contains(b)
                }
                AlignmentCoordinates::Secondary { descendant, .. } => {
                    (..end.secondary_ordinate_descendant().unwrap()).contains(descendant)
                }
            }
        } else {
            true
        }
    }

    pub fn can_increment_both(&self, end: AlignmentCoordinates) -> bool {
        self.can_increment_a(end) && self.can_increment_b(end)
    }

    pub fn increment_a(&self) -> Self {
        match self {
            AlignmentCoordinates::Primary { a, b } => {
                AlignmentCoordinates::Primary { a: a + 1, b: *b }
            }
            AlignmentCoordinates::Secondary {
                ancestor,
                descendant,
                ts_kind,
            } => AlignmentCoordinates::Secondary {
                ancestor: ancestor - 1,
                descendant: *descendant,
                ts_kind: *ts_kind,
            },
        }
    }

    pub fn increment_b(&self) -> Self {
        match self {
            AlignmentCoordinates::Primary { a, b } => {
                AlignmentCoordinates::Primary { a: *a, b: b + 1 }
            }
            AlignmentCoordinates::Secondary {
                ancestor,
                descendant,
                ts_kind,
            } => AlignmentCoordinates::Secondary {
                ancestor: *ancestor,
                descendant: descendant + 1,
                ts_kind: *ts_kind,
            },
        }
    }

    pub fn increment_both(&self) -> Self {
        self.increment_a().increment_b()
    }

    /// Generate all possible 12-jumps.
    ///
    /// The TS kind is given by the start coordinates.
    /// The left and right limits of the jump are given by the start and end coordinates.
    /// The end coordinates must be in primary form and simply be the end of the aligned sequences.
    /// The start coordinates are in secondary form.
    pub fn generate_12_jumps(
        &self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> impl Iterator<Item = (isize, Self)> {
        let Self::Primary { a, b } = *self else {
            panic!("Can only generate 12-jumps from primary coordinates");
        };
        let ts_kind = start.ts_kind().unwrap();
        let ancestor_zero = match ts_kind.ancestor {
            TsAncestor::Seq1 => a,
            TsAncestor::Seq2 => b,
        } as isize;
        let ancestor_limit = match ts_kind.ancestor {
            TsAncestor::Seq1 => end.primary_ordinate_a().unwrap(),
            TsAncestor::Seq2 => end.primary_ordinate_b().unwrap(),
        };
        let descendant = match ts_kind.descendant {
            TsDescendant::Seq1 => a,
            TsDescendant::Seq2 => b,
        };

        (start.secondary_ordinate_ancestor().unwrap()..ancestor_limit).map(move |ancestor| {
            (
                ancestor as isize - ancestor_zero,
                Self::Secondary {
                    ancestor,
                    descendant,
                    ts_kind,
                },
            )
        })
    }
}

impl Display for AlignmentCoordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentCoordinates::Primary { a, b } => write!(f, "({}, {})", a, b),
            AlignmentCoordinates::Secondary {
                ancestor,
                descendant,
                ts_kind,
            } => write!(f, "({}, {}, {:?})", ancestor, descendant, ts_kind),
        }
    }
}
