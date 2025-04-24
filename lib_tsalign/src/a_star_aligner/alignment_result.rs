use std::fmt::{Display, Formatter, Result, Write};

use a_star_sequences::SequencePair;
use alignment::Alignment;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStarResult, cost::AStarCost};
use noisy_float::types::{R64, r64};
use num_traits::{Float, Zero};

pub mod a_star_sequences;
pub mod alignment;

pub trait IAlignmentType {
    fn is_repeatable(&self) -> bool;

    fn is_repeated(&self, previous: &Self) -> bool;

    fn is_internal(&self) -> bool;

    fn is_template_switch_exit(&self) -> bool;
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(tag = "type"))]
pub enum AlignmentResult<AlignmentType, Cost> {
    WithTarget {
        #[cfg_attr(feature = "serde", serde(flatten))]
        alignment: Alignment<AlignmentType>,

        #[cfg_attr(feature = "serde", serde(flatten))]
        statistics: AlignmentStatistics<Cost>,
    },

    WithoutTarget {
        #[cfg_attr(feature = "serde", serde(flatten))]
        statistics: AlignmentStatistics<Cost>,
    },
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[must_use]
pub struct AlignmentStatistics<Cost> {
    pub result: AStarResult<(), Cost>,

    // Input data used by show.
    pub sequences: SequencePair,
    pub reference_offset: usize,
    pub query_offset: usize,

    pub cost: R64,
    pub cost_per_base: R64,
    pub duration_seconds: R64,
    pub opened_nodes: R64,
    pub closed_nodes: R64,
    pub suboptimal_opened_nodes: R64,
    pub suboptimal_opened_nodes_ratio: R64,
    pub template_switch_amount: R64,

    /// Runtime in seconds.
    ///
    /// To be filled by some other tool, not collected by tsalign.
    #[cfg_attr(feature = "serde", serde(default))]
    pub runtime: R64,

    /// Memory in bytes.
    ///
    /// To be filled by some other tool, not collected by tsalign.
    #[cfg_attr(feature = "serde", serde(default))]
    pub memory: R64,
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
        $action!(template_switch_amount);
        $action!(runtime);
        $action!(memory);
    }};
}

impl<AlignmentType: IAlignmentType, Cost: AStarCost> AlignmentResult<AlignmentType, Cost> {
    #[expect(clippy::too_many_arguments)]
    pub fn new_with_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        alignment: Vec<(usize, AlignmentType)>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        reference_offset: usize,
        query_offset: usize,
        result: AStarResult<(), Cost>,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self::new(
            Some(alignment),
            reference,
            query,
            reference_name,
            query_name,
            reference_offset,
            query_offset,
            result,
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            reference_length,
            query_length,
        )
    }

    #[expect(clippy::too_many_arguments)]
    pub fn new_without_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        result: AStarResult<(), Cost>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        reference_offset: usize,
        query_offset: usize,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self::new(
            None,
            reference,
            query,
            reference_name,
            query_name,
            reference_offset,
            query_offset,
            result,
            duration_seconds,
            opened_nodes,
            closed_nodes,
            suboptimal_opened_nodes,
            reference_length,
            query_length,
        )
    }

    #[expect(clippy::too_many_arguments)]
    fn new<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        alignment: Option<Vec<(usize, AlignmentType)>>,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        reference_name: &str,
        query_name: &str,
        reference_offset: usize,
        query_offset: usize,
        result: AStarResult<(), Cost>,
        duration_seconds: f64,
        opened_nodes: usize,
        closed_nodes: usize,
        suboptimal_opened_nodes: usize,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        let cost = result.cost();
        let statistics = AlignmentStatistics {
            result,
            sequences: SequencePair::new(reference, query, reference_name, query_name),
            reference_offset,
            query_offset,

            cost: (cost.as_f64()).try_into().unwrap(),
            cost_per_base: ((cost.as_f64() * 2.0) / (reference_length + query_length) as f64)
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
            template_switch_amount: r64(alignment
                .as_ref()
                .map(|alignment| {
                    alignment
                        .iter()
                        .filter(|(_, alignment_type)| alignment_type.is_template_switch_exit())
                        .count() as f64
                })
                .unwrap_or(0.0)),
            runtime: r64(0.0),
            memory: r64(0.0),
        };

        if let Some(alignment) = alignment {
            Self::WithTarget {
                alignment: alignment.into(),
                statistics,
            }
        } else {
            Self::WithoutTarget { statistics }
        }
    }
}

impl<AlignmentType, Cost> AlignmentResult<AlignmentType, Cost> {
    pub fn statistics(&self) -> &AlignmentStatistics<Cost> {
        match self {
            AlignmentResult::WithTarget { statistics, .. } => statistics,
            AlignmentResult::WithoutTarget { statistics } => statistics,
        }
    }

    pub fn statistics_mut(&mut self) -> &mut AlignmentStatistics<Cost> {
        match self {
            AlignmentResult::WithTarget { statistics, .. } => statistics,
            AlignmentResult::WithoutTarget { statistics } => statistics,
        }
    }
}

impl<AlignmentType: IAlignmentType, Cost> AlignmentResult<AlignmentType, Cost> {
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
        let Self::WithTarget { alignment, .. } = self else {
            return Ok(());
        };

        alignment.write_cigar(writer)
    }
}

impl<Cost: Clone> AlignmentStatistics<Cost> {
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
        assert!(!statistics.is_empty());
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

impl<AlignmentType: Display + IAlignmentType, Cost: Display> Display
    for AlignmentResult<AlignmentType, Cost>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if let AlignmentResult::WithTarget { .. } = self {
            write!(f, "CIGAR: ")?;
            self.write_cigar(f)?;
            writeln!(f)?;
        } else {
            writeln!(f, "No alignment found")?;
        }

        let (AlignmentResult::WithTarget { statistics, .. }
        | AlignmentResult::WithoutTarget { statistics }) = self;

        statistics.fmt(f)
    }
}

impl<Cost: Display> Display for AlignmentStatistics<Cost> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        writeln!(f, "{}", self.result)?;
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

impl<Cost> Default for AlignmentStatistics<Cost> {
    fn default() -> Self {
        Self {
            result: Default::default(),
            sequences: Default::default(),
            reference_offset: Default::default(),
            query_offset: Default::default(),
            cost: Default::default(),
            cost_per_base: Default::default(),
            duration_seconds: Default::default(),
            opened_nodes: Default::default(),
            closed_nodes: Default::default(),
            suboptimal_opened_nodes: Default::default(),
            suboptimal_opened_nodes_ratio: Default::default(),
            template_switch_amount: Default::default(),
            runtime: Default::default(),
            memory: Default::default(),
        }
    }
}
