use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{Context, NodeData};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait NodeOrdStrategy: AlignmentStrategy {
    fn cmp(&self, n1: &NodeData, n2: &NodeData) -> std::cmp::Ordering;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct CostOnlyNodeOrdStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AntiDiagonalNodeOrdStrategy;

impl NodeOrdStrategy for CostOnlyNodeOrdStrategy {
    fn cmp(&self, n1: &NodeData, n2: &NodeData) -> std::cmp::Ordering {
        n1.lower_bound_cost().cmp(&n2.lower_bound_cost())
    }
}

impl NodeOrdStrategy for AntiDiagonalNodeOrdStrategy {
    fn cmp(&self, n1: &NodeData, n2: &NodeData) -> std::cmp::Ordering {
        match n1.lower_bound_cost().cmp(&n2.lower_bound_cost()) {
            // This secondary ordering may make things actually slower.
            // While it does reduce the number of visited nodes a little bit,
            // it also makes heap operations more expensive.
            // Preliminary testing showed that this would be a slowdown.
            std::cmp::Ordering::Equal => n2
                .identifier
                .anti_diagonal()
                .cmp(&n1.identifier.anti_diagonal()),
            ordering => ordering,
        }
    }
}

impl AlignmentStrategy for CostOnlyNodeOrdStrategy {
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

impl AlignmentStrategy for AntiDiagonalNodeOrdStrategy {
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
