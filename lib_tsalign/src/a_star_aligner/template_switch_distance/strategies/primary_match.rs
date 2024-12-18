use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::cost::Cost;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier};

use super::AlignmentStrategySelector;

pub trait PrimaryMatchStrategy: Clone + Eq + Debug + Display {
    type Memory;
    type Identifier: Copy + Debug + Ord + Hash;

    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        identifier: Identifier<Strategies>,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn create_root_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier;

    fn generate_successor_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier;

    /// If false, then this node cannot have a primary match outgoing edge outside of a flank.
    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    /// If false, then this node cannot have a primary match outgoing edge within a flank.
    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool;

    fn fake_substitution_cost<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        &self,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Cost;

    fn available_primary_matches(&self, identifier: Self::Identifier) -> usize;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AllowPrimaryMatchStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MaxConsecutivePrimaryMatchStrategy;

pub struct MaxConsecutivePrimaryMatchMemory {
    pub max_consecutive_primary_matches: usize,
    pub root_available_primary_matches: usize,
    pub fake_substitution_cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct MaxConsecutivePrimaryMatchIdentifier {
    available_primary_matches: usize,
}

impl PrimaryMatchStrategy for AllowPrimaryMatchStrategy {
    type Memory = ();
    type Identifier = ();

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
        _identifier: Identifier<Strategies>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }

    fn create_root_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier {
    }

    fn generate_successor_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier {
        identifier
    }

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Self::Identifier,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        // Can always use a match.
        true
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Self::Identifier,
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

    fn available_primary_matches(&self, _identifier: Self::Identifier) -> usize {
        usize::MAX
    }
}

impl PrimaryMatchStrategy for MaxConsecutivePrimaryMatchStrategy {
    type Memory = MaxConsecutivePrimaryMatchMemory;
    type Identifier = MaxConsecutivePrimaryMatchIdentifier;

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
        _identifier: Identifier<Strategies>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }

    fn create_root_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier {
        Self::Identifier {
            available_primary_matches: context.memory.primary_match.root_available_primary_matches,
        }
    }

    fn generate_successor_identifier<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self::Identifier {
        Self::Identifier {
            available_primary_matches: if alignment_type == AlignmentType::PrimaryMatch {
                debug_assert!(identifier.available_primary_matches > 0);
                identifier.available_primary_matches - 1
            } else {
                context.memory.primary_match.max_consecutive_primary_matches
            },
        }
    }

    fn can_do_primary_non_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        identifier: Self::Identifier,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> bool {
        debug_assert!(identifier.available_primary_matches < isize::MAX as usize);
        identifier.available_primary_matches > 0
    }

    fn can_do_primary_flank_match<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<PrimaryMatch = Self>,
    >(
        _identifier: Self::Identifier,
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

    fn available_primary_matches(&self, identifier: Self::Identifier) -> usize {
        identifier.available_primary_matches
    }
}

impl Display for AllowPrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrimaryMatch allowed")
    }
}

impl Display for MaxConsecutivePrimaryMatchStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrimaryMatch")
    }
}
