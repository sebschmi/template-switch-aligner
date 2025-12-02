use std::fmt::Display;

use crate::alignment::TsKind;

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

    pub fn can_increment_a(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        assert_eq!(self.ts_kind(), start.ts_kind());
        assert_eq!(self.ts_kind(), end.ts_kind());

        match self {
            AlignmentCoordinates::Primary { a, .. } => {
                (start.primary_ordinate_a().unwrap()..end.primary_ordinate_a().unwrap()).contains(a)
            }
            AlignmentCoordinates::Secondary { ancestor, .. } => {
                (start.secondary_ordinate_ancestor().unwrap()
                    ..end.secondary_ordinate_ancestor().unwrap())
                    .contains(ancestor)
            }
        }
    }

    pub fn can_increment_b(&self, start: AlignmentCoordinates, end: AlignmentCoordinates) -> bool {
        assert_eq!(self.ts_kind(), start.ts_kind());
        assert_eq!(self.ts_kind(), end.ts_kind());

        match self {
            AlignmentCoordinates::Primary { b, .. } => {
                (start.primary_ordinate_b().unwrap()..end.primary_ordinate_b().unwrap()).contains(b)
            }
            AlignmentCoordinates::Secondary { descendant, .. } => {
                (start.secondary_ordinate_descendant().unwrap()
                    ..end.secondary_ordinate_descendant().unwrap())
                    .contains(descendant)
            }
        }
    }

    pub fn can_increment_both(
        &self,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> bool {
        self.can_increment_a(start, end) && self.can_increment_b(start, end)
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
