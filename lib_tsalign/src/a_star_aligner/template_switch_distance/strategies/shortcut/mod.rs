use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::AStarContext;

use crate::a_star_aligner::template_switch_distance::{
    lower_bounds::template_switch::TemplateSwitchLowerBoundMatrix, Context,
};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait ShortcutStrategy: AlignmentStrategy {
    type Memory;

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Shortcut = Self>,
    >(
        node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        context: &Context<SubsequenceType, Strategies>,
        opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    );
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoShortcutStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct TemplateSwitchShortcutStrategy;

impl ShortcutStrategy for NoShortcutStrategy {
    type Memory = ();

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Shortcut = Self>,
    >(
        _node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        _context: &Context<SubsequenceType, Strategies>,
        _opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    ) {
        // Do nothing.
    }
}

impl ShortcutStrategy for TemplateSwitchShortcutStrategy {
    type Memory = TemplateSwitchLowerBoundMatrix;

    fn generate_successors<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector<Shortcut = Self>,
    >(
        node: &<Context<SubsequenceType, Strategies> as AStarContext>::Node,
        context: &Context<SubsequenceType, Strategies>,
        opened_nodes_output: &mut impl for<'reference, 'query> Extend<
            <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
        >,
    ) {
        todo!()
    }
}

impl AlignmentStrategy for NoShortcutStrategy {
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

impl AlignmentStrategy for TemplateSwitchShortcutStrategy {
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
