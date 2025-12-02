//! Representation of an alignment.

use std::fmt::Display;

pub mod coordinates;
pub mod sequences;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum AlignmentType {
    Match,
    Substitution,
    GapA,
    GapB,
    TsStart { jump: isize, ts_kind: TsKind },
    TsEnd,
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

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct TsKind {
    ancestor: TsAncestor,
    descendant: TsDescendant,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum GapType {
    None,
    /// A gap in sequence 1, meaning that sequence 2 has characters that are missing from sequence 1.
    InA,
    /// A gap in sequence 2, meaning that sequence 1 has characters that are missing from sequence 2.
    InB,
}

pub struct Alignment {
    pub alignment: Vec<(usize, AlignmentType)>,
}

impl Display for GapType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GapType::None => write!(f, "M/S"),
            GapType::InA => write!(f, "GA"),
            GapType::InB => write!(f, "GB"),
        }
    }
}

impl FromIterator<AlignmentType> for Alignment {
    fn from_iter<T: IntoIterator<Item = AlignmentType>>(iter: T) -> Self {
        let mut alignment = Vec::new();
        for alignment_type in iter {
            if Some(alignment_type) == alignment.last().map(|(_, alignment_type)| *alignment_type) {
                alignment.last_mut().unwrap().0 += 1;
            } else {
                alignment.push((1, alignment_type));
            }
        }
        Self { alignment }
    }
}

impl From<Vec<AlignmentType>> for Alignment {
    fn from(alignment_types: Vec<AlignmentType>) -> Self {
        alignment_types.into_iter().collect()
    }
}
