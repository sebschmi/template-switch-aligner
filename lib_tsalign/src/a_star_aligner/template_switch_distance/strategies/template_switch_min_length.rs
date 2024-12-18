use std::{collections::HashMap, marker::PhantomData, mem};

use compact_genome::interface::sequence::GenomeSequence;
use deterministic_default_hasher::DeterministicDefaultHasher;
use generic_a_star::reset::Reset;
use generic_a_star::{AStar, AStarContext, AStarNode, AStarResult};

use crate::a_star_aligner::template_switch_distance::AlignmentType;
use crate::{
    a_star_aligner::template_switch_distance::{
        identifier::{GapType, TemplateSwitchPrimary, TemplateSwitchSecondary},
        Context, Identifier, Node,
    },
    costs::cost::Cost,
};

use super::primary_match::PrimaryMatchStrategy;
use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait TemplateSwitchMinLengthStrategy: AlignmentStrategy {
    /// The type used to memorise lookahead results.
    type Memory: Default + Reset;

    /// Takes the template switch entrance node and provides a lower bound for its costs depending on the minimum length of a template switch.
    /// The modified entrance node is returned in the iterator along with further nodes that were created while computing the lower bound.
    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector<TemplateSwitchMinLength = Self>,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        secondary_root_node: Node<Strategies>,
        context: &mut Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>>;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoTemplateSwitchMinLengthStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct LookaheadTemplateSwitchMinLengthStrategy;

struct TemplateSwitchMinLengthContext<
    'reference,
    'query,
    'context,
    SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    Strategies: AlignmentStrategySelector,
> {
    context: &'context mut Context<'reference, 'query, SubsequenceType, Strategies>,
    root_node: Node<Strategies>,
    phantom_data: PhantomData<Strategies>,
}

impl TemplateSwitchMinLengthStrategy for NoTemplateSwitchMinLengthStrategy {
    type Memory = ();

    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector<TemplateSwitchMinLength = Self>,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        secondary_root_node: Node<Strategies>,
        _context: &mut Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        Some(secondary_root_node)
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct LookaheadMemoryKey {
    template_switch_primary: TemplateSwitchPrimary,
    template_switch_secondary: TemplateSwitchSecondary,
    primary_index: usize,
    secondary_index: usize,
}

impl TemplateSwitchMinLengthStrategy for LookaheadTemplateSwitchMinLengthStrategy {
    type Memory = HashMap<LookaheadMemoryKey, Cost, DeterministicDefaultHasher>;

    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector<TemplateSwitchMinLength = Self>,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        mut secondary_root_node: Node<Strategies>,
        context: &mut Context<SubsequenceType, Strategies>,
    ) -> impl IntoIterator<Item = Node<Strategies>> {
        let Identifier::Secondary {
            template_switch_primary,
            template_switch_secondary,
            length: 0,
            primary_index,
            secondary_index,
            gap_type: GapType::None,
            ..
        } = secondary_root_node.node_data.identifier
        else {
            unreachable!("Only called with a secondary root node.")
        };

        let memory_key = LookaheadMemoryKey {
            template_switch_primary,
            template_switch_secondary,
            primary_index,
            secondary_index,
        };

        let secondary_root_node = if let Some(a_star_lower_bound) =
            context.memory.template_switch_min_length.get(&memory_key)
        {
            secondary_root_node.node_data.a_star_lower_bound += *a_star_lower_bound;
            secondary_root_node
        } else {
            let buffers = mem::take(&mut context.a_star_buffers);
            let initial_cost = secondary_root_node.cost();
            let mut a_star = AStar::new_with_buffers(
                TemplateSwitchMinLengthContext::new(secondary_root_node.clone(), context),
                buffers,
            );
            a_star.initialise();

            let alignment_result = a_star.search();

            if let AStarResult::FoundTarget { identifier, .. } = alignment_result {
                let target_cost = a_star.closed_node(&identifier).unwrap().cost();
                context.a_star_buffers = a_star.into_buffers();
                let lower_bound = target_cost - initial_cost;

                context
                    .memory
                    .template_switch_min_length
                    .insert(memory_key, lower_bound);
                secondary_root_node.node_data.a_star_lower_bound += lower_bound;

                secondary_root_node
            } else {
                context.a_star_buffers = a_star.into_buffers();
                secondary_root_node
            }
        };

        Some(secondary_root_node)
    }
}

impl AlignmentStrategy for NoTemplateSwitchMinLengthStrategy {
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
        _identifier: Identifier<<<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy>::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl AlignmentStrategy for LookaheadTemplateSwitchMinLengthStrategy {
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
        _identifier: Identifier<<<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy>::IdentifierPrimaryExtraData>,
        _alignment_type: AlignmentType,
        _context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        *self
    }
}

impl<
        'reference,
        'query,
        'context,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > TemplateSwitchMinLengthContext<'reference, 'query, 'context, SubsequenceType, Strategies>
{
    fn new(
        root_node: Node<Strategies>,
        context: &'context mut Context<'reference, 'query, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            context,
            root_node,
            phantom_data: PhantomData,
        }
    }
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > AStarContext for TemplateSwitchMinLengthContext<'_, '_, '_, SubsequenceType, Strategies>
{
    type Node = Node<Strategies>;

    fn create_root(&self) -> Self::Node {
        self.root_node.clone()
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        self.context.generate_successors(node, output);
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        let Identifier::Secondary { length, .. } = node.node_data.identifier else {
            unreachable!("A non-secondary node was closed before a target was closed.")
        };
        debug_assert!(length <= self.context.config.min_length);

        length == self.context.config.min_length
    }
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > Reset for TemplateSwitchMinLengthContext<'_, '_, '_, SubsequenceType, Strategies>
{
    fn reset(&mut self) {
        unimplemented!("Designed to be used only once.")
    }
}