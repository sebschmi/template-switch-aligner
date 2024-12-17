use std::fmt::{Debug, Display};

use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::cost::Cost;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::AlignmentStrategySelector;

pub trait PrimaryMatchStrategy: Eq + Clone + Debug + Display {
    type Memory;

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _identifier: Identifier,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

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

    fn available_primary_matches(&self) -> usize;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AllowPrimaryMatchStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxConsecutivePrimaryMatchStrategy {
    available_primary_matches: usize,
}

pub struct MaxConsecutivePrimaryMatchMemory {
    pub max_consecutive_primary_matches: usize,
    pub root_available_primary_matches: usize,
    pub fake_substitution_cost: Cost,
}

impl PrimaryMatchStrategy for AllowPrimaryMatchStrategy {
    type Memory = ();

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _identifier: Identifier,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }

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

    fn available_primary_matches(&self) -> usize {
        usize::MAX
    }
}

impl PrimaryMatchStrategy for MaxConsecutivePrimaryMatchStrategy {
    type Memory = MaxConsecutivePrimaryMatchMemory;

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            available_primary_matches: context.memory.primary_match.root_available_primary_matches,
        }
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _identifier: Identifier,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            available_primary_matches: if alignment_type == AlignmentType::PrimaryMatch {
                debug_assert!(self.available_primary_matches > 0);
                self.available_primary_matches - 1
            } else {
                context.memory.primary_match.max_consecutive_primary_matches
            },
        }
    }

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        debug_assert!(self.available_primary_matches < isize::MAX as usize);
        self.available_primary_matches > 0
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

    fn available_primary_matches(&self) -> usize {
        self.available_primary_matches
    }
}

impl Display for AllowPrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrimaryMatch allowed")
    }
}

impl Display for MaxConsecutivePrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PrimaryMatch available: {}",
            self.available_primary_matches
        )
    }
}
