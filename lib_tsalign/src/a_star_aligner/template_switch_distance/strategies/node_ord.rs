use compact_genome::interface::sequence::GenomeSequence;

use crate::a_star_aligner::template_switch_distance::{AlignmentType, Context, Identifier, Node};

use super::{
    primary_match::{MaxConsecutivePrimaryMatchStrategy, PrimaryMatchStrategy},
    AlignmentStrategy, AlignmentStrategySelector,
};

pub trait NodeOrdStrategy<PrimaryMatch: PrimaryMatchStrategy>: AlignmentStrategy {
    fn cmp<Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct CostOnlyNodeOrdStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct AntiDiagonalNodeOrdStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct PrimaryMatchNodeOrdStrategy;

impl<PrimaryMatch: PrimaryMatchStrategy> NodeOrdStrategy<PrimaryMatch> for CostOnlyNodeOrdStrategy {
    fn cmp<Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering {
        n1.node_data
            .lower_bound_cost()
            .cmp(&n2.node_data.lower_bound_cost())
    }
}

impl<PrimaryMatch: PrimaryMatchStrategy> NodeOrdStrategy<PrimaryMatch>
    for AntiDiagonalNodeOrdStrategy
{
    fn cmp<Strategies: AlignmentStrategySelector<PrimaryMatch = PrimaryMatch>>(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering {
        match n1
            .node_data
            .lower_bound_cost()
            .cmp(&n2.node_data.lower_bound_cost())
        {
            // This secondary ordering may make things actually slower.
            // While it does reduce the number of visited nodes a little bit,
            // it also makes heap operations more expensive.
            // Preliminary testing showed that this would be a slowdown.
            std::cmp::Ordering::Equal => n2
                .node_data
                .identifier
                .anti_diagonal()
                .cmp(&n1.node_data.identifier.anti_diagonal()),
            ordering => ordering,
        }
    }
}

impl NodeOrdStrategy<MaxConsecutivePrimaryMatchStrategy> for PrimaryMatchNodeOrdStrategy {
    fn cmp<
        Strategies: AlignmentStrategySelector<PrimaryMatch = MaxConsecutivePrimaryMatchStrategy>,
    >(
        &self,
        n1: &Node<Strategies>,
        n2: &Node<Strategies>,
    ) -> std::cmp::Ordering {
        n1.node_data
            .lower_bound_cost()
            .cmp(&n2.node_data.lower_bound_cost())
            .then_with(|| {
                n2.strategies
                    .primary_match
                    .available_primary_matches()
                    .cmp(&n1.strategies.primary_match.available_primary_matches())
            })
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
        _identifier: Identifier<()>,
        _alignment_type: AlignmentType,
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
        _identifier: Identifier<()>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
impl AlignmentStrategy for PrimaryMatchNodeOrdStrategy {
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
        _identifier: Identifier<()>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}
