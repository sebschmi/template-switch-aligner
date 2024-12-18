use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::cost::Cost;
use log::debug;
use num_traits::SaturatingSub;
use seed_chain::{
    chain::{Chain, ChainingCostsProvider},
    seed::{ChainingAnchor, ChainingAnchors},
};

use crate::{
    a_star_aligner::template_switch_distance::{
        identifier::GapType,
        lower_bounds::{
            template_switch::TemplateSwitchLowerBoundMatrix,
            template_switch_alignment::TemplateSwitchAlignmentLowerBoundMatrix,
        },
        AlignmentType, Context, Identifier, Node,
    },
    config::TemplateSwitchConfig,
};

use super::{primary_match::PrimaryMatchStrategy, AlignmentStrategy, AlignmentStrategySelector};

pub trait ChainingStrategy: AlignmentStrategy {
    type Memory;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType>,
        block_size: usize,
    ) -> Self::Memory;

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        nodes: impl IntoIterator<Item = Node<Strategies>>,
        context: &Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>>;
}

#[expect(dead_code)]
pub struct ChainingMemory {
    ts_lower_bounds: TemplateSwitchLowerBoundMatrix,
    tsa_lower_bounds: TemplateSwitchAlignmentLowerBoundMatrix,
    chain: Chain,
    max_gap_open_cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoChainingStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PrecomputeOnlyChainingStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct LowerBoundChainingStrategy;

struct TemplateSwitchAlignmentLowerBoundChainingCosts<'a> {
    matrix: &'a TemplateSwitchAlignmentLowerBoundMatrix,
    reference_length: usize,
    query_length: usize,
}

impl ChainingStrategy for NoChainingStrategy {
    type Memory = ();

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        _reference: &SubsequenceType,
        _query: &SubsequenceType,
        _config: &TemplateSwitchConfig<AlphabetType>,
        _block_size: usize,
    ) -> Self::Memory {
        // Do nothing.
    }

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        nodes: impl IntoIterator<Item = Node<Strategies>>,
        _context: &Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        nodes
    }
}

impl ChainingStrategy for PrecomputeOnlyChainingStrategy {
    type Memory = ChainingMemory;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType>,
        block_size: usize,
    ) -> Self::Memory {
        let ts_lower_bounds = TemplateSwitchLowerBoundMatrix::new(config);
        debug!("{ts_lower_bounds}");
        let tsa_lower_bounds = TemplateSwitchAlignmentLowerBoundMatrix::new(
            config,
            &ts_lower_bounds,
            reference.len(),
            query.len(),
            block_size * 2 - 1,
            block_size - 1,
        );
        debug!("{tsa_lower_bounds}");
        let chaining_anchors = ChainingAnchors::seed_nonoverlapping(reference, query, block_size);
        let chain = Chain::compute_chain(
            TemplateSwitchAlignmentLowerBoundChainingCosts {
                matrix: &tsa_lower_bounds,
                reference_length: reference.len(),
                query_length: query.len(),
            },
            chaining_anchors,
        );
        debug!("{chain}");

        ChainingMemory {
            ts_lower_bounds,
            tsa_lower_bounds,
            chain,
            max_gap_open_cost: config.primary_edit_costs.max_gap_open_cost(),
        }
    }

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        nodes: impl IntoIterator<Item = Node<Strategies>>,
        _context: &Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        nodes
    }
}

impl ChainingStrategy for LowerBoundChainingStrategy {
    type Memory = ChainingMemory;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType>,
        block_size: usize,
    ) -> Self::Memory {
        PrecomputeOnlyChainingStrategy::initialise_memory(reference, query, config, block_size)
    }

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        nodes: impl IntoIterator<Item = Node<Strategies>>,
        context: &Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        nodes.into_iter().map(|mut node| {
            if let Identifier::Primary {
                reference_index,
                query_index,
                gap_type,
                flank_index,
                ..
            }
            | Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type,
                flank_index,
                ..
            } = node.node_data.identifier
            {
                if flank_index <= 0 {
                    let mut chain_lower_bound = context
                        .memory
                        .chaining
                        .chain
                        .chain_lower_bound(reference_index, query_index);
                    if gap_type != GapType::None {
                        chain_lower_bound = chain_lower_bound
                            .saturating_sub(&context.memory.chaining.max_gap_open_cost);
                    }

                    node.node_data.a_star_lower_bound =
                        node.node_data.a_star_lower_bound.max(chain_lower_bound);
                }
            }

            node
        })
    }
}

impl ChainingCostsProvider for TemplateSwitchAlignmentLowerBoundChainingCosts<'_> {
    fn chaining_costs(
        &self,
        from: &seed_chain::chain::Identifier,
        to: &seed_chain::chain::Identifier,
    ) -> Cost {
        let from = seed_chain_identifier_to_chaining_anchor(
            from,
            self.reference_length,
            self.query_length,
        );
        let to =
            seed_chain_identifier_to_chaining_anchor(to, self.reference_length, self.query_length);
        if from.reference_block().end > to.reference_block().start
            || from.query_block().end > to.query_block().start
        {
            return Cost::MAX;
        }

        let delta_reference = to.reference_block().start - from.reference_block().end;
        let delta_query = to.query_block().start - from.query_block().end;
        self.matrix.cost(delta_reference, delta_query)
    }
}

fn seed_chain_identifier_to_chaining_anchor(
    identifier: &seed_chain::chain::Identifier,
    reference_length: usize,
    query_length: usize,
) -> ChainingAnchor {
    match identifier {
        seed_chain::chain::Identifier::Root => ChainingAnchor::new(0..0, 0..0),
        seed_chain::chain::Identifier::Anchor { anchor } => anchor.clone(),
        seed_chain::chain::Identifier::Target => ChainingAnchor::new(
            reference_length..reference_length,
            query_length..query_length,
        ),
    }
}

impl AlignmentStrategy for NoChainingStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<<<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy>::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl AlignmentStrategy for PrecomputeOnlyChainingStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<<<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy>::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl AlignmentStrategy for LowerBoundChainingStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<<<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy>::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
