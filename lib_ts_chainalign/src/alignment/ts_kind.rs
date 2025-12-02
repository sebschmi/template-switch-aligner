use std::fmt::Display;

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct TsKind {
    pub ancestor: TsAncestor,
    pub descendant: TsDescendant,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum TsAncestor {
    Seq1,
    Seq2,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum TsDescendant {
    Seq1,
    Seq2,
}

impl TsKind {
    pub const TS11: Self = TsKind {
        ancestor: TsAncestor::Seq1,
        descendant: TsDescendant::Seq1,
    };
    pub const TS12: Self = TsKind {
        ancestor: TsAncestor::Seq1,
        descendant: TsDescendant::Seq2,
    };
    pub const TS21: Self = TsKind {
        ancestor: TsAncestor::Seq2,
        descendant: TsDescendant::Seq1,
    };
    pub const TS22: Self = TsKind {
        ancestor: TsAncestor::Seq2,
        descendant: TsDescendant::Seq2,
    };
}

impl Display for TsKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TS{}{}",
            match self.ancestor {
                TsAncestor::Seq1 => "1",
                TsAncestor::Seq2 => "2",
            },
            match self.descendant {
                TsDescendant::Seq1 => "1",
                TsDescendant::Seq2 => "2",
            }
        )
    }
}
