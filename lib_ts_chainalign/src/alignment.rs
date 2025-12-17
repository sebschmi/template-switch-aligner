//! Representation of an alignment.

use std::{fmt::Display, iter};

use crate::alignment::ts_kind::TsKind;

pub mod coordinates;
pub mod sequences;
pub mod ts_kind;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum AlignmentType {
    Match,
    Substitution,
    GapA,
    GapB,
    TsStart { jump: isize, ts_kind: TsKind },
    TsEnd { jump: isize },
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum GapType {
    None,
    /// A gap in sequence 1, meaning that sequence 2 has characters that are missing from sequence 1.
    InA,
    /// A gap in sequence 2, meaning that sequence 1 has characters that are missing from sequence 2.
    InB,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub alignment: Vec<(usize, AlignmentType)>,
}

impl Alignment {
    pub fn iter_unpacked(&self) -> impl Iterator<Item = AlignmentType> {
        self.alignment
            .iter()
            .flat_map(|(multiplicity, alignment_type)| {
                iter::repeat_n(*alignment_type, *multiplicity)
            })
    }
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

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (multiplicity, alignment_type) in &self.alignment {
            write!(f, "{multiplicity}{alignment_type}")?;
        }
        Ok(())
    }
}

impl Display for AlignmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentType::Match => write!(f, "M"),
            AlignmentType::Substitution => write!(f, "S"),
            AlignmentType::GapA => write!(f, "GA"),
            AlignmentType::GapB => write!(f, "GB"),
            AlignmentType::TsStart { jump, ts_kind } => write!(f, "TS{}[{jump}]", ts_kind.digits()),
            AlignmentType::TsEnd { jump } => write!(f, "END[{jump}]"),
        }
    }
}
