use std::fmt::Debug;

use compact_genome::{
    implementation::{
        alphabets::{
            dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
            dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet,
            rna_alphabet::RnaAlphabet, rna_alphabet_or_n::RnaAlphabetOrN,
            rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet,
        },
        vec_sequence::VectorGenome,
    },
    interface::{
        alphabet::Alphabet,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
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
            node_ord::{AntiDiagonalNodeOrdStrategy, NodeOrdStrategy},
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

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub struct Config {
    pub alphabet: InputAlphabet,
    pub reference_name: String,
    pub query_name: String,
    /// Costs specification in plain text format.
    ///
    /// This is the same format that the `.tsa` config files use.
    pub costs: String,

    pub node_ord_strategy: NodeOrdStrategySelector,
    pub min_length_strategy: MinLengthStrategySelector,
    pub chaining_strategy: ChainingStrategySelector,
    pub total_length_strategy: TotalLengthStrategySelector,
    pub no_ts: bool,

    pub cost_limit: Option<U64Cost>,
    /// Approximate memory limit in bytes.
    pub memory_limit: Option<usize>,
    pub range: Option<AlignmentRange>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            alphabet: InputAlphabet::DnaN,
            reference_name: "reference".to_owned(),
            query_name: "query".to_owned(),
            costs: TemplateSwitchConfig::<DnaAlphabetOrN, U64Cost>::default().to_string(),
            node_ord_strategy: NodeOrdStrategySelector::AntiDiagonal,
            min_length_strategy: MinLengthStrategySelector::Lookahead,
            chaining_strategy: ChainingStrategySelector::None,
            total_length_strategy: TotalLengthStrategySelector::Maximise,
            no_ts: false,
            cost_limit: None,
            memory_limit: None,
            range: None,
        }
    }
}

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum InputAlphabet {
    Dna,
    DnaN,
    Rna,
    RnaN,
    DnaIupac,
    RnaIupac,
}

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum NodeOrdStrategySelector {
    CostOnly,
    AntiDiagonal,
}

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum MinLengthStrategySelector {
    None,
    Lookahead,
}

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum ChainingStrategySelector {
    None,
    PrecomputeOnly,
    LowerBound,
}

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "snake_case"))]
pub enum TotalLengthStrategySelector {
    None,
    Maximise,
}

/// Align `query` to `reference` with the given `config`.
///
/// `query` and `reference` must be ASCII strings restricted to the characters specified by `config.alphabet`.
pub fn a_star_align(
    reference: &[u8],
    query: &[u8],
    config: &Config,
) -> AlignmentResult<AlignmentType, U64Cost> {
    match config.alphabet {
        InputAlphabet::Dna => {
            a_star_align_select_node_ord_strategy::<DnaAlphabet>(reference, query, config)
        }
        InputAlphabet::DnaN => {
            a_star_align_select_node_ord_strategy::<DnaAlphabetOrN>(reference, query, config)
        }
        InputAlphabet::Rna => {
            a_star_align_select_node_ord_strategy::<RnaAlphabet>(reference, query, config)
        }
        InputAlphabet::RnaN => {
            a_star_align_select_node_ord_strategy::<RnaAlphabetOrN>(reference, query, config)
        }
        InputAlphabet::DnaIupac => a_star_align_select_node_ord_strategy::<
            DnaIupacNucleicAcidAlphabet,
        >(reference, query, config),
        InputAlphabet::RnaIupac => a_star_align_select_node_ord_strategy::<
            RnaIupacNucleicAcidAlphabet,
        >(reference, query, config),
    }
}

fn a_star_align_select_node_ord_strategy<AlphabetType: Alphabet + Debug + Clone + Eq>(
    reference: &[u8],
    query: &[u8],
    config: &Config,
) -> AlignmentResult<AlignmentType, U64Cost> {
    let reference = VectorGenome::<AlphabetType>::from_slice_u8(reference).unwrap();
    let query = VectorGenome::from_slice_u8(query).unwrap();
    let costs = TemplateSwitchConfig::read_plain(config.costs.as_bytes()).unwrap();

    match config.node_ord_strategy {
        NodeOrdStrategySelector::CostOnly => {
            /*a_star_align_select_template_switch_min_length_strategy::<_, _, CostOnlyNodeOrdStrategy>(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                config,
                costs,
            )*/
            unimplemented!("The other option appears to always be better.");
        }
        NodeOrdStrategySelector::AntiDiagonal => {
            a_star_align_select_template_switch_min_length_strategy::<
                _,
                _,
                AntiDiagonalNodeOrdStrategy,
            >(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                config,
                costs,
            )
        }
    }
}

fn a_star_align_select_template_switch_min_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: &Config,
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,
) -> AlignmentResult<AlignmentType, U64Cost> {
    match config.min_length_strategy {
        MinLengthStrategySelector::None => a_star_align_select_chaining_strategy::<
            _,
            _,
            NodeOrd,
            NoTemplateSwitchMinLengthStrategy<U64Cost>,
        >(reference, query, config, costs),
        MinLengthStrategySelector::Lookahead => a_star_align_select_chaining_strategy::<
            _,
            _,
            NodeOrd,
            LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
        >(reference, query, config, costs),
    }
}

fn a_star_align_select_chaining_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: &Config,
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,
) -> AlignmentResult<AlignmentType, U64Cost> {
    match config.chaining_strategy {
        ChainingStrategySelector::None => a_star_align_select_no_ts_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            NoChainingStrategy<U64Cost>,
        >(reference, query, config, costs),
        ChainingStrategySelector::PrecomputeOnly => {
            /*a_star_align_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy<U64Cost>,
            >(reference, query, config, costs)*/
            unimplemented!("No reason to precompute without using the information.");
        }
        ChainingStrategySelector::LowerBound => a_star_align_select_no_ts_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            LowerBoundChainingStrategy<U64Cost>,
        >(reference, query, config, costs),
    }
}

fn a_star_align_select_no_ts_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: &Config,
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,
) -> AlignmentResult<AlignmentType, U64Cost> {
    if config.no_ts {
        a_star_align_select_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            MaxTemplateSwitchCountStrategy,
        >(reference, query, config, costs, 0)
    } else {
        a_star_align_select_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
        >(reference, query, config, costs, ())
    }
}

fn a_star_align_select_total_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: &Config,
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) -> AlignmentResult<AlignmentType, U64Cost> {
    match config.total_length_strategy {
        TotalLengthStrategySelector::None => a_star_align_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            TemplateSwitchCount,
            NoTemplateSwitchTotalLengthStrategy,
        >(
            reference,
            query,
            config,
            costs,
            template_switch_count_memory,
        ),
        TotalLengthStrategySelector::Maximise => a_star_align_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            TemplateSwitchCount,
            MaxTemplateSwitchTotalLengthStrategy,
        >(
            reference,
            query,
            config,
            costs,
            template_switch_count_memory,
        ),
    }
}

fn a_star_align_call<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    TemplateSwitchTotalLength: TemplateSwitchTotalLengthStrategy,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: &Config,
    costs: TemplateSwitchConfig<AlphabetType, U64Cost>,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) -> AlignmentResult<AlignmentType, U64Cost> {
    template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            AlphabetType,
            U64Cost,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            TemplateSwitchCount,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            TemplateSwitchTotalLength,
        >,
        _,
    >(
        reference,
        query,
        &config.reference_name,
        &config.query_name,
        config.range.clone(),
        costs,
        config.cost_limit,
        config.memory_limit,
        false,
        template_switch_count_memory,
    )
}
