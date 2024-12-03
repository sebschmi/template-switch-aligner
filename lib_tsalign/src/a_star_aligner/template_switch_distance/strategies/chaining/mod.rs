use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use log::debug;

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

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory;
}

#[expect(dead_code)]
pub struct ChainingMemory {
    ts_lower_bounds: TemplateSwitchLowerBoundMatrix,
    tsa_lower_bounds: TemplateSwitchAlignmentLowerBoundMatrix,
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
        // Do nothing.
    }
}

impl ChainingStrategy for PrecomputeOnlyChainingStrategy {
    type Memory = ChainingMemory;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory {
        let ts_lower_bounds = TemplateSwitchLowerBoundMatrix::new(config);
        debug!("{ts_lower_bounds}");
        let tsa_lower_bounds =
            TemplateSwitchAlignmentLowerBoundMatrix::new(config, &ts_lower_bounds, 100, 100, 10);
        debug!("{tsa_lower_bounds}");

        ChainingMemory {
            ts_lower_bounds,
            tsa_lower_bounds,
        }
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
