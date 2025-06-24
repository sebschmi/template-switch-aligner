use std::{fmt::Debug, marker::PhantomData};

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::cost::AStarCost;
use get_size2::GetSize;
use log::debug;
use seed_chain::{
    chain::{Chain, ChainingCostsProvider},
    seed::{ChainingAnchor, ChainingAnchors},
};

use crate::{
    a_star_aligner::template_switch_distance::{
        AlignmentType, Context, Identifier, Node,
        identifier::GapType,
        lower_bounds::{
            template_switch::TemplateSwitchLowerBoundMatrix,
            template_switch_alignment::TemplateSwitchAlignmentLowerBoundMatrix,
        },
    },
    config::TemplateSwitchConfig,
};

use super::{AlignmentStrategy, AlignmentStrategySelector, primary_match::PrimaryMatchStrategy};

pub trait ChainingStrategy<Cost>: AlignmentStrategy {
    type Memory: GetSize;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
        block_size: usize,
    ) -> Self::Memory;

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Cost = Cost, Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        node: Node<Strategies>,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Node<Strategies>;
}

#[derive(GetSize)]
#[allow(dead_code)] // The compiler won't recognize this as dead code since the memory estimation is using it, but it factually is dead code.
pub struct ChainingMemory<Cost> {
    ts_lower_bounds: TemplateSwitchLowerBoundMatrix<Cost>,
    tsa_lower_bounds: TemplateSwitchAlignmentLowerBoundMatrix<Cost>,
    chain: Chain<Cost>,
    max_gap_open_cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoChainingStrategy<Cost> {
    phantom_data: PhantomData<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PrecomputeOnlyChainingStrategy<Cost> {
    phantom_data: PhantomData<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct LowerBoundChainingStrategy<Cost> {
    phantom_data: PhantomData<Cost>,
}

#[derive(GetSize)]
struct TemplateSwitchAlignmentLowerBoundChainingCosts<'a, Cost> {
    matrix: &'a TemplateSwitchAlignmentLowerBoundMatrix<Cost>,
    reference_length: usize,
    query_length: usize,
    phantom_data: PhantomData<Cost>,
}

impl<Cost: AStarCost> ChainingStrategy<Cost> for NoChainingStrategy<Cost> {
    type Memory = ();

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        _reference: &SubsequenceType,
        _query: &SubsequenceType,
        _config: &TemplateSwitchConfig<AlphabetType, Cost>,
        _block_size: usize,
    ) -> Self::Memory {
        // Do nothing.
    }

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Cost = Cost, Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        node: Node<Strategies>,
        _context: &Context<SubsequenceType, Strategies>,
    ) -> Node<Strategies> {
        node
    }
}

impl<Cost: AStarCost> ChainingStrategy<Cost> for PrecomputeOnlyChainingStrategy<Cost> {
    type Memory = ChainingMemory<Cost>;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
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
                phantom_data: PhantomData,
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
        Strategies: AlignmentStrategySelector<Cost = Cost, Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        node: Node<Strategies>,
        _context: &Context<SubsequenceType, Strategies>,
    ) -> Node<Strategies> {
        node
    }
}

impl<Cost: AStarCost> ChainingStrategy<Cost> for LowerBoundChainingStrategy<Cost> {
    type Memory = ChainingMemory<Cost>;

    fn initialise_memory<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        config: &TemplateSwitchConfig<AlphabetType, Cost>,
        block_size: usize,
    ) -> Self::Memory {
        PrecomputeOnlyChainingStrategy::initialise_memory(reference, query, config, block_size)
    }

    fn apply_lower_bound<
        Strategies: AlignmentStrategySelector<Cost = Cost, Chaining = Self>,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        mut node: Node<Strategies>,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Node<Strategies> {
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
    }
}

impl<Cost: AStarCost> ChainingCostsProvider
    for TemplateSwitchAlignmentLowerBoundChainingCosts<'_, Cost>
{
    type Cost = Cost;

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
            return Cost::max_value();
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

impl<Cost: AStarCost> AlignmentStrategy for NoChainingStrategy<Cost> {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            phantom_data: PhantomData,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl<Cost: AStarCost> AlignmentStrategy for PrecomputeOnlyChainingStrategy<Cost> {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            phantom_data: PhantomData,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl<Cost: AStarCost> AlignmentStrategy for LowerBoundChainingStrategy<Cost> {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            phantom_data: PhantomData,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
