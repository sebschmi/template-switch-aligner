use std::fmt::{Display, Formatter, Result, Write};

use noisy_float::types::R64;
use num_traits::{Bounded, Zero};

use crate::costs::cost::Cost;

pub trait IAlignmentType {
    fn is_repeatable(&self) -> bool;

    fn is_repeated(&self, previous: &Self) -> bool;

    fn is_internal(&self) -> bool;
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AlignmentResult<AlignmentType> {
    pub alignment: Vec<(usize, AlignmentType)>,
    pub statistics: AlignmentStatistics,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AlignmentStatistics {
    pub cost: R64,
    pub cost_per_base: R64,
    pub duration_seconds: R64,
    pub opened_nodes: R64,
    pub closed_nodes: R64,
    pub suboptimal_opened_nodes: R64,
    pub suboptimal_opened_nodes_ratio: R64,
}

impl<AlignmentType> AlignmentResult<AlignmentType> {
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        alignment: Vec<(usize, AlignmentType)>,
        cost: Cost,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self {
            alignment,
            statistics: AlignmentStatistics {
                cost: (cost.as_u64() as f64).try_into().unwrap(),
                cost_per_base: ((cost.as_u64() * 2) as f64
                    / (reference_length + query_length) as f64)
                    .try_into()
                    .unwrap(),
                duration_seconds: duration_seconds.try_into().unwrap(),
                opened_nodes: (opened_nodes as f64).try_into().unwrap(),
                closed_nodes: (closed_nodes as f64).try_into().unwrap(),
                suboptimal_opened_nodes: (suboptimal_opened_nodes as f64).try_into().unwrap(),
                suboptimal_opened_nodes_ratio: (suboptimal_opened_nodes as f64
                    / (opened_nodes - suboptimal_opened_nodes) as f64)
                    .try_into()
                    .unwrap(),
            },
        }
    }
}

impl<AlignmentType: IAlignmentType> AlignmentResult<AlignmentType> {
    pub fn cigar(&self) -> String
    where
        AlignmentType: Display,
    {
        let mut result = String::new();
        self.write_cigar(&mut result).unwrap();
        result
    }

    pub fn write_cigar(&self, writer: &mut impl Write) -> Result
    where
        AlignmentType: Display,
    {
        for (amount, alignment_type) in &self.alignment {
            if alignment_type.is_repeatable() {
                write!(writer, "{amount}{alignment_type}")?;
            } else {
                write!(writer, "{alignment_type}")?;
            }
        }

        Ok(())
    }
}

impl AlignmentStatistics {
    pub fn min_value() -> Self {
        AlignmentStatistics {
            cost: R64::min_value(),
            cost_per_base: R64::min_value(),
            duration_seconds: R64::min_value(),
            opened_nodes: R64::min_value(),
            closed_nodes: R64::min_value(),
            suboptimal_opened_nodes: R64::min_value(),
            suboptimal_opened_nodes_ratio: R64::min_value(),
        }
    }

    pub fn max_value() -> Self {
        AlignmentStatistics {
            cost: R64::max_value(),
            cost_per_base: R64::max_value(),
            duration_seconds: R64::max_value(),
            opened_nodes: R64::max_value(),
            closed_nodes: R64::max_value(),
            suboptimal_opened_nodes: R64::max_value(),
            suboptimal_opened_nodes_ratio: R64::max_value(),
        }
    }

    pub fn zero_value() -> Self {
        AlignmentStatistics {
            cost: R64::zero(),
            cost_per_base: R64::zero(),
            duration_seconds: R64::zero(),
            opened_nodes: R64::zero(),
            closed_nodes: R64::zero(),
            suboptimal_opened_nodes: R64::zero(),
            suboptimal_opened_nodes_ratio: R64::zero(),
        }
    }
}

impl<AlignmentType: Display + IAlignmentType> Display for AlignmentResult<AlignmentType> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        write!(f, "CIGAR: ")?;
        self.write_cigar(f)?;
        writeln!(f)?;
        self.statistics.fmt(f)
    }
}

impl Display for AlignmentStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        writeln!(f, "Cost: {}", self.cost)?;
        writeln!(f, "Cost per base: {:.2}", self.cost_per_base)?;
        writeln!(f, "Opened nodes: {}", self.opened_nodes)?;
        writeln!(f, "Closed nodes: {}", self.closed_nodes)?;
        writeln!(
            f,
            "Suboptimal openend nodes: {}",
            self.suboptimal_opened_nodes
        )?;
        writeln!(
            f,
            "Suboptimal openend nodes per optimal opened node: {:.2}",
            self.suboptimal_opened_nodes_ratio
        )?;
        write!(f, "Duration: {:.2}s", self.duration_seconds)?;

        Ok(())
    }
}
