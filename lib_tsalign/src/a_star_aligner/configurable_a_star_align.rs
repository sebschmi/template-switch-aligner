use std::fmt::Debug;

use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
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
            chaining::{
                ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy,
                PrecomputeOnlyChainingStrategy,
            },
            node_ord::{AntiDiagonalNodeOrdStrategy, CostOnlyNodeOrdStrategy, NodeOrdStrategy},
            primary_match::AllowPrimaryMatchStrategy,
            template_switch_count::{
                MaxTemplateSwitchCountStrategy, NoTemplateSwitchCountStrategy,
                TemplateSwitchCountStrategy,
            },
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                TemplateSwitchMinLengthStrategy,
            },
        },
    },
};

#[derive(Debug)]
pub struct Config {
    alphabet: InputAlphabet,
    reference_name: String,
    query_name: String,
    /// Costs specification in plain text format.
    ///
    /// This is the same format that the `.tsa` config files use.
    costs: String,

    node_ord_strategy: NodeOrdStrategySelector,
    min_length_strategy: MinLengthStrategySelector,
    chaining_strategy: ChainingStrategySelector,
    no_ts: bool,

    cost_limit: Option<U64Cost>,
    /// Approximate memory limit in bytes.
    memory_limit: Option<usize>,
    range: Option<AlignmentRange>,
}

#[derive(Debug)]
pub enum InputAlphabet {
    Dna,
    DnaN,
    Rna,
    RnaN,
    DnaIupac,
    RnaIupac,
}

#[derive(Debug)]
pub enum NodeOrdStrategySelector {
    CostOnly,
    AntiDiagonal,
}

#[derive(Debug)]
pub enum MinLengthStrategySelector {
    None,
    Lookahead,
}

#[derive(Debug)]
pub enum ChainingStrategySelector {
    None,
    PrecomputeOnly,
    LowerBound,
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
        InputAlphabet::DnaN => todo!(),
        InputAlphabet::Rna => todo!(),
        InputAlphabet::RnaN => todo!(),
        InputAlphabet::DnaIupac => todo!(),
        InputAlphabet::RnaIupac => todo!(),
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
            a_star_align_select_template_switch_min_length_strategy::<_, _, CostOnlyNodeOrdStrategy>(
                reference.as_genome_subsequence(),
                query.as_genome_subsequence(),
                config,
                costs,
            )
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
        MinLengthStrategySelector::None => align_a_star_template_switch_select_chaining_strategy::<
            _,
            _,
            NodeOrd,
            NoTemplateSwitchMinLengthStrategy<U64Cost>,
        >(reference, query, config, costs),
        MinLengthStrategySelector::Lookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            >(reference, query, config, costs)
        }
    }
}

fn align_a_star_template_switch_select_chaining_strategy<
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
        ChainingStrategySelector::None => align_a_star_template_switch_select_no_ts_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            NoChainingStrategy<U64Cost>,
        >(reference, query, config, costs),
        ChainingStrategySelector::PrecomputeOnly => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy<U64Cost>,
            >(reference, query, config, costs)
        }
        ChainingStrategySelector::LowerBound => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                LowerBoundChainingStrategy<U64Cost>,
            >(reference, query, config, costs)
        }
    }
}

fn align_a_star_template_switch_select_no_ts_strategy<
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
        align_a_star_template_switch_distance_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            MaxTemplateSwitchCountStrategy,
        >(reference, query, config, costs, 0)
    } else {
        align_a_star_template_switch_distance_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
        >(reference, query, config, costs, ())
    }
}

fn align_a_star_template_switch_distance_call<
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
        template_switch_count_memory,
    )
}
