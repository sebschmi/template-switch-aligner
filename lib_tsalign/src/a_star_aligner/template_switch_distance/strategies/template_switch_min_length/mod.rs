use crate::a_star_aligner::template_switch_distance::{Context, Node};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait TemplateSwitchMinLengthStrategy: AlignmentStrategy {
    /// Takes the template switch entrance node and provides a lower bound for its costs depending on the minimum length of a template switch.
    /// The modified entrance node is returned in the iterator along with further nodes that were created while computing the lower bound.
    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        secondary_root_node: Node<Strategies>,
        context: &mut Context<Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>>;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoTemplateSwitchMinLengthStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct LookaheadTemplateSwitchMinLengthStrategy;

impl TemplateSwitchMinLengthStrategy for NoTemplateSwitchMinLengthStrategy {
    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        _reference: &SubsequenceType,
        _query: &SubsequenceType,
        secondary_root_node: Node<Strategies>,
        _context: &mut Context<Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        Some(secondary_root_node)
    }
}

impl TemplateSwitchMinLengthStrategy for LookaheadTemplateSwitchMinLengthStrategy {
    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        secondary_root_node: Node<Strategies>,
        context: &mut Context<Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        todo!();
        []
    }
}

impl AlignmentStrategy for NoTemplateSwitchMinLengthStrategy {
    fn create_root<Strategies: AlignmentStrategySelector>(_context: &Context<Strategies>) -> Self {
        Self
    }

    fn generate_successor<Strategies: AlignmentStrategySelector>(
        &self,
        _context: &Context<Strategies>,
    ) -> Self {
        *self
    }
}

impl AlignmentStrategy for LookaheadTemplateSwitchMinLengthStrategy {
    fn create_root<Strategies: AlignmentStrategySelector>(_context: &Context<Strategies>) -> Self {
        Self
    }

    fn generate_successor<Strategies: AlignmentStrategySelector>(
        &self,
        _context: &Context<Strategies>,
    ) -> Self {
        *self
    }
}
