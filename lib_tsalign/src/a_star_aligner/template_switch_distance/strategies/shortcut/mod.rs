use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStarContext, AStarNode};

use crate::{
    a_star_aligner::template_switch_distance::{
        lower_bounds::template_switch::TemplateSwitchLowerBoundMatrix, Context, Identifier,
    },
    config::TemplateSwitchConfig,
};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait ShortcutStrategy: AlignmentStrategy {
    type Memory;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory;

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
pub struct TemplateSwitchLowerBoundShortcutStrategy;

impl ShortcutStrategy for NoShortcutStrategy {
    type Memory = ();

    fn initialise_memory<AlphabetType: Alphabet>(
        _config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory {
        // Do nothing.
    }

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

impl ShortcutStrategy for TemplateSwitchLowerBoundShortcutStrategy {
    type Memory = TemplateSwitchLowerBoundMatrix;

    fn initialise_memory<AlphabetType: Alphabet>(
        config: &TemplateSwitchConfig<AlphabetType>,
    ) -> Self::Memory {
        TemplateSwitchLowerBoundMatrix::new(config)
    }

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
        match *node.identifier() {
            Identifier::Primary { flank_index, .. }
            | Identifier::PrimaryReentry { flank_index, .. } => {
                if flank_index == context.config.left_flank_length {
                    opened_nodes_output.extend(context.memory.shortcut.iter().flat_map(|entry| {
                        node.generate_primary_shortcut_successor(
                            entry.x(),
                            entry.y(),
                            entry.cost(),
                            context,
                        )
                    }));
                }
            }

            _ => { /* Do nothing. */ }
        }
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

impl AlignmentStrategy for TemplateSwitchLowerBoundShortcutStrategy {
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
