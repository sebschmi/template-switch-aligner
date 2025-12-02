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
