use std::{fmt::Debug, str::FromStr};

use compact_genome::{
    implementation::{alphabets::dna_alphabet_or_n::DnaAlphabetOrN, vec_sequence::VectorGenome},
    interface::sequence::{GenomeSequence, OwnedGenomeSequence},
};
use generic_a_star::cost::U64Cost;

use crate::{
    a_star_aligner::{
        template_switch_distance::strategies::{
            AlignmentStrategySelection, primary_range::NoPrunePrimaryRangeStrategy,
            secondary_deletion::AllowSecondaryDeletionStrategy, shortcut::NoShortcutStrategy,
        },
        template_switch_distance_a_star_align,
    },
    config::TemplateSwitchConfig,
};

use super::{
    alignment_geometry::AlignmentRange,
    alignment_result::AlignmentResult,
    template_switch_distance::{
        AlignmentType,
        strategies::{
            chaining::{ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy},
            node_ord::AntiDiagonalNodeOrdStrategy,
            primary_match::AllowPrimaryMatchStrategy,
            template_switch_count::{
                MaxTemplateSwitchCountStrategy, NoTemplateSwitchCountStrategy,
                TemplateSwitchCountStrategy,
            },
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                TemplateSwitchMinLengthStrategy,
            },
            template_switch_total_length::{
                MaxTemplateSwitchTotalLengthStrategy, NoTemplateSwitchTotalLengthStrategy,
                TemplateSwitchTotalLengthStrategy,
            },
        },
    },
};

pub use compact_genome::implementation::alphabets;
pub use compact_genome::interface::alphabet::Alphabet;

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum MinLengthStrategySelector {
    None,
    #[default]
    Lookahead,
}

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum ChainingStrategySelector {
    #[default]
    None,
    // PrecomputeOnly,
    LowerBound,
}

#[derive(Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum TotalLengthStrategySelector {
    None,
    #[default]
    Maximise,
}

/// Just used in this file to bundle the query parameters to make the code more readable
struct QueryData<'a> {
    reference_name: &'a str,
    reference: &'a [u8],
    query_name: &'a str,
    query: &'a [u8],
    ranges: Option<AlignmentRange>,
    cost_limit: Option<u64>,
    memory_limit: Option<usize>,
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize))] // Mostly used for kwargs in python bindings
#[cfg_attr(feature = "serde", serde(default))]
pub struct Aligner<AlphabetType: Alphabet = DnaAlphabetOrN> {
    #[cfg_attr(feature = "serde", serde(skip))] // Not deserializable
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,

    // ↓ Settings for how the alignment (definition)
    min_length_strategy: MinLengthStrategySelector,
    chaining_strategy: ChainingStrategySelector,
    total_length_strategy: TotalLengthStrategySelector,
    no_ts: bool,
}

impl Aligner<DnaAlphabetOrN> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// We implement this ourselves to avoid the overly restrictive trait bound
/// `AlphabetType: Default` that would be generated with the derive macro
impl<AlphabetType: Alphabet> Default for Aligner<AlphabetType> {
    fn default() -> Self {
        Self {
            costs: TemplateSwitchConfig::<AlphabetType, U64Cost>::default(),
            min_length_strategy: MinLengthStrategySelector::default(),
            chaining_strategy: ChainingStrategySelector::default(),
            total_length_strategy: TotalLengthStrategySelector::default(),
            no_ts: false,
        }
    }
}

impl<AlphabetType: Alphabet> Aligner<AlphabetType> {
    /// Parse a cost string and update the aligner cost
    ///
    /// # Errors
    /// Returns an error if the string cannot be parsed as `TemplateSwitchConfig<AlphabetType>`
    pub fn set_costs_parse(&mut self, costs: &str) -> Result<&mut Self, crate::error::Error> {
        self.costs = TemplateSwitchConfig::from_str(costs)?;
        Ok(self)
    }

    pub fn set_costs(&mut self, costs: TemplateSwitchConfig<AlphabetType, U64Cost>) -> &mut Self {
        self.costs = costs;
        self
    }

    pub fn set_min_length_strategy(
        &mut self,
        min_length_strategy: MinLengthStrategySelector,
    ) -> &mut Self {
        self.min_length_strategy = min_length_strategy;
        self
    }

    pub fn set_chaining_strategy(
        &mut self,
        chaining_strategy: ChainingStrategySelector,
    ) -> &mut Self {
        self.chaining_strategy = chaining_strategy;
        self
    }

    pub fn set_total_length_strategy(
        &mut self,
        total_length_strategy: TotalLengthStrategySelector,
    ) -> &mut Self {
        self.total_length_strategy = total_length_strategy;
        self
    }

    pub fn set_no_ts(&mut self, no_ts: bool) -> &mut Self {
        self.no_ts = no_ts;
        self
    }

    /// This performs the actual alignment and is, depending on the inputs, potentially taking quite long.
    /// Consider spawning tasks to not block the user application.
    #[allow(
        clippy::too_many_arguments,
        reason = "For now, this is the intended API"
    )]
    #[must_use]
    pub fn align(
        &self,
        reference_name: &str,
        reference: &[u8],
        query_name: &str,
        query: &[u8],
        ranges: Option<AlignmentRange>,
        cost_limit: Option<u64>,
        memory_limit: Option<usize>,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        let data = QueryData {
            reference_name,
            reference,
            query_name,
            query,
            ranges,
            cost_limit,
            memory_limit,
        };
        self.align_select_min_length_strategy(data)
    }

    fn align_select_min_length_strategy(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.min_length_strategy {
            MinLengthStrategySelector::None => self
                .align_select_chaining_strategy::<NoTemplateSwitchMinLengthStrategy<U64Cost>>(data),
            MinLengthStrategySelector::Lookahead => self.align_select_chaining_strategy::<LookaheadTemplateSwitchMinLengthStrategy<U64Cost>>(data),
        }
    }

    fn align_select_chaining_strategy<ML: TemplateSwitchMinLengthStrategy<U64Cost>>(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.chaining_strategy {
            ChainingStrategySelector::None => {
                self.align_select_total_length_strategy::<ML, NoChainingStrategy<U64Cost>>(data)
            }
            ChainingStrategySelector::LowerBound => self
                .align_select_total_length_strategy::<ML, LowerBoundChainingStrategy<U64Cost>>(
                    data,
                ),
        }
    }

    fn align_select_total_length_strategy<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
    >(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        match self.total_length_strategy {
            TotalLengthStrategySelector::None => self
                .align_select_no_ts_strategy::<ML, CH, NoTemplateSwitchTotalLengthStrategy>(data),
            TotalLengthStrategySelector::Maximise => self
                .align_select_no_ts_strategy::<ML, CH, MaxTemplateSwitchTotalLengthStrategy>(data),
        }
    }

    fn align_select_no_ts_strategy<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
    >(
        &self,
        data: QueryData,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        if self.no_ts {
            self.align_call::<ML, CH, TL, MaxTemplateSwitchCountStrategy>(data, 0)
        } else {
            self.align_call::<ML, CH, TL, NoTemplateSwitchCountStrategy>(data, ())
        }
    }

    fn align_call<
        ML: TemplateSwitchMinLengthStrategy<U64Cost>,
        CH: ChainingStrategy<U64Cost>,
        TL: TemplateSwitchTotalLengthStrategy,
        TC: TemplateSwitchCountStrategy,
    >(
        &self,
        data: QueryData,
        count_strategy_memory: TC::Memory,
    ) -> AlignmentResult<AlignmentType, U64Cost> {
        // TODO error handling
        let reference = VectorGenome::<AlphabetType>::from_slice_u8(data.reference).unwrap();
        let query = VectorGenome::from_slice_u8(data.query).unwrap();
        let cost_limit = data.cost_limit.map(U64Cost::from);
        template_switch_distance_a_star_align::<
            AlignmentStrategySelection<
                AlphabetType,
                U64Cost,
                AntiDiagonalNodeOrdStrategy,
                ML,
                CH,
                TC,
                AllowSecondaryDeletionStrategy,
                NoShortcutStrategy<U64Cost>,
                AllowPrimaryMatchStrategy,
                NoPrunePrimaryRangeStrategy,
                TL,
            >,
            _,
        >(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            data.reference_name,
            data.query_name,
            data.ranges,
            &self.costs,
            cost_limit,
            data.memory_limit,
            false,
            count_strategy_memory,
        )
    }
}
