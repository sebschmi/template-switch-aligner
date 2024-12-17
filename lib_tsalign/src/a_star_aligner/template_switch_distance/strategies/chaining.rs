use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::cost::Cost;
use log::debug;
use seed_chain::{
    chain::{Chain, ChainingCostsProvider},
    seed::{ChainingAnchor, ChainingAnchors},
};

use crate::{
    a_star_aligner::template_switch_distance::{
        lower_bounds::{
            template_switch::TemplateSwitchLowerBoundMatrix,
            template_switch_alignment::TemplateSwitchAlignmentLowerBoundMatrix,
        },
        AlignmentType, Context, Identifier,
    },
    config::TemplateSwitchConfig,
};

use super::{AlignmentStrategy, AlignmentStrategySelector};

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
}

#[expect(dead_code)]
pub struct ChainingMemory {
    ts_lower_bounds: TemplateSwitchLowerBoundMatrix,
    tsa_lower_bounds: TemplateSwitchAlignmentLowerBoundMatrix,
    chain: Chain,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoChainingStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PrecomputeOnlyChainingStrategy;

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
        }
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
        _identifier: Identifier,
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
        _identifier: Identifier,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
