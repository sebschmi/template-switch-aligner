use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use log::info;

use crate::{
    a_star_aligner::template_switch_distance::Context, config::TemplateSwitchConfig,
    template_switch_lower_bound::TemplateSwitchLowerBoundMatrix,
};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait ChainingStrategy: AlignmentStrategy {
    type Memory;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory;
}

#[expect(dead_code)]
pub struct ChainingMemory {
    ts_lower_bounds: TemplateSwitchLowerBoundMatrix,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoChainingStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PrecomputeOnlyChainingStrategy;

impl ChainingStrategy for NoChainingStrategy {
    type Memory = ();

    fn initialise_memory<AlphabetType: Alphabet>(
        _config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory {
    }
}

impl ChainingStrategy for PrecomputeOnlyChainingStrategy {
    type Memory = ChainingMemory;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory {
        info!("Computing TS lower bound matrix...");
        let ts_lower_bounds = TemplateSwitchLowerBoundMatrix::new(config);
        println!("{ts_lower_bounds}");
        ChainingMemory { ts_lower_bounds }
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
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
