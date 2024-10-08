use std::fmt::{Display, Formatter, Result, Write};

use noisy_float::types::R64;
use num_traits::{Float, Zero};

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

    #[cfg_attr(feature = "serde", serde(flatten))]
    pub statistics: AlignmentStatistics,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[must_use]
pub struct AlignmentStatistics {
    pub cost: R64,
    pub cost_per_base: R64,
    pub duration_seconds: R64,
    pub opened_nodes: R64,
    pub closed_nodes: R64,
    pub suboptimal_opened_nodes: R64,
    pub suboptimal_opened_nodes_ratio: R64,
}

macro_rules! each_statistic {
    ($action:path) => {{
        $action!(cost);
        $action!(cost_per_base);
        $action!(duration_seconds);
        $action!(opened_nodes);
        $action!(closed_nodes);
        $action!(suboptimal_opened_nodes);
        $action!(suboptimal_opened_nodes_ratio);
    }};
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
        let mut result = Self::default();

        macro_rules! create_min_value {
            ($field:ident) => {
                result.$field = R64::min_value();
            };
        }
        each_statistic!(create_min_value);

        result
    }

    pub fn max_value() -> Self {
        let mut result = Self::default();

        macro_rules! create_max_value {
            ($field:ident) => {
                result.$field = R64::max_value();
            };
        }
        each_statistic!(create_max_value);

        result
    }

    pub fn zero() -> Self {
        let mut result = Self::default();

        macro_rules! create_zero {
            ($field:ident) => {
                result.$field = R64::zero();
            };
        }
        each_statistic!(create_zero);

        result
    }

    pub fn piecewise_min(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! piecewise_min {
            ($field:ident) => {
                result.$field = result.$field.min(other.$field);
            };
        }
        each_statistic!(piecewise_min);

        result
    }

    pub fn piecewise_max(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! piecewise_max {
            ($field:ident) => {
                result.$field = result.$field.max(other.$field);
            };
        }
        each_statistic!(piecewise_max);

        result
    }

    pub fn piecewise_add(&self, other: &Self) -> Self {
        let mut result = self.clone();

        macro_rules! add {
            ($field:ident) => {
                result.$field += other.$field;
            };
        }
        each_statistic!(add);

        result
    }

    pub fn piecewise_div(&self, divisor: R64) -> Self {
        let mut result = self.clone();

        macro_rules! div {
            ($field:ident) => {
                result.$field /= divisor;
            };
        }
        each_statistic!(div);

        result
    }

    pub fn piecewise_percentile(statistics: &[Self], percentile: R64) -> Self {
        assert!(percentile >= 0.0);
        assert!(percentile <= 1.0);
        let mut result = Self::zero();
        let mut buffer = Vec::new();

        macro_rules! percentile {
            ($field:ident) => {
                buffer.clear();
                buffer.extend(statistics.iter().map(|s| s.$field));
                buffer.sort();

                // Scale percentile.
                let index = (percentile * buffer.len() as f64).floor().raw() as usize;
                let index = if index == buffer.len() {
                    // Edge case if percentile == 1.0
                    index - 1
                } else {
                    index
                };
                result.$field = buffer[index];
            };
        }
        each_statistic!(percentile);

        result
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
