use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::cost::Cost;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait PrimaryMatchStrategy: AlignmentStrategy {
    type Memory;

    /// If false, then this node cannot have a primary match outgoing edge outside of a flank.
    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    /// If false, then this node cannot have a primary match outgoing edge within a flank.
    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AllowPrimaryMatchStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxConsecutivePrimaryMatchStrategy {
    consecutive_primary_matches: usize,
}

pub struct MaxConsecutivePrimaryMatchMemory {
    pub max_consecutive_primary_matches: usize,
    pub fake_substitution_cost: Cost,
}

impl PrimaryMatchStrategy for AllowPrimaryMatchStrategy {
    type Memory = ();

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a match.
        true
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a match.
        true
    }

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost {
        Cost::MAX
    }
}

impl PrimaryMatchStrategy for MaxConsecutivePrimaryMatchStrategy {
    type Memory = MaxConsecutivePrimaryMatchMemory;

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        debug_assert!(
            self.consecutive_primary_matches
                <= context.memory.primary_match.max_consecutive_primary_matches
        );
        self.consecutive_primary_matches
            < context.memory.primary_match.max_consecutive_primary_matches
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a flank match.
        true
    }

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost {
        context.memory.primary_match.fake_substitution_cost
    }
}

impl AlignmentStrategy for AllowPrimaryMatchStrategy {
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

impl AlignmentStrategy for MaxConsecutivePrimaryMatchStrategy {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            consecutive_primary_matches: 0,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        _identifier: Identifier,
        alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            consecutive_primary_matches: if alignment_type == AlignmentType::PrimaryMatch {
                self.consecutive_primary_matches + 1
            } else {
                0
            },
        }
    }
}
