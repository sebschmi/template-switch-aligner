//! Representation of an alignment.

use std::fmt::Display;

pub mod coordinates;
pub mod sequences;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum AlignmentType {
    Match,
    Substitution,
    Gap1,
    Gap2,
    TsStart,
    TsEnd,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum GapType {
    None,
    /// A gap in sequence 1, meaning that sequence 2 has characters that are missing from sequence 1.
    In1,
    /// A gap in sequence 2, meaning that sequence 1 has characters that are missing from sequence 2.
    In2,
}

pub struct Alignment {
    pub alignment: Vec<(usize, AlignmentType)>,
}

impl Display for GapType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GapType::None => write!(f, "M/S"),
            GapType::In1 => write!(f, "GA"),
            GapType::In2 => write!(f, "GB"),
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
