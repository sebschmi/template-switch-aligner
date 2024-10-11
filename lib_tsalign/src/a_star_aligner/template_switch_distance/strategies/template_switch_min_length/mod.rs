use std::{collections::HashMap, mem};

use binary_heap_plus::BinaryHeap;

use crate::{
    a_star_aligner::{
        a_star_align_loop,
        template_switch_distance::{
            identifier::{GapType, TemplateSwitchPrimary, TemplateSwitchSecondary},
            Context, Identifier, Node,
        },
        AlignmentGraphNode,
    },
    costs::cost::Cost,
};

use super::{AlignmentStrategy, AlignmentStrategySelector};

pub trait TemplateSwitchMinLengthStrategy: AlignmentStrategy {
    /// The type used to memorise lookahead results.
    type Memory: Default;

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
        reference: &SubsequenceType,
        query: &SubsequenceType,
        secondary_root_node: Node<Strategies>,
        context: &mut Context<Strategies>,
        closed_nodes_output: &mut impl Extend<(
            <Node<Strategies> as AlignmentGraphNode<Strategies::Alphabet>>::Identifier,
            Node<Strategies>,
        )>,
    ) -> impl IntoIterator<Item = Node<Strategies>>;
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct NoTemplateSwitchMinLengthStrategy;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct LookaheadTemplateSwitchMinLengthStrategy;

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
        _reference: &SubsequenceType,
        _query: &SubsequenceType,
        secondary_root_node: Node<Strategies>,
        _context: &mut Context<Strategies>,
        _closed_nodes_output: &mut impl Extend<(
            <Node<Strategies> as AlignmentGraphNode<Strategies::Alphabet>>::Identifier,
            Node<Strategies>,
        )>,
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
    type Memory = HashMap<LookaheadMemoryKey, Cost>;

    fn template_switch_min_length_lookahead<
        Strategies: AlignmentStrategySelector<TemplateSwitchMinLength = Self>,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        mut secondary_root_node: Node<Strategies>,
        context: &mut Context<Strategies>,
        closed_nodes_output: &mut impl Extend<(
            <Node<Strategies> as AlignmentGraphNode<Strategies::Alphabet>>::Identifier,
            Node<Strategies>,
        )>,
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

        if let Some(a_star_lower_bound) = context.template_switch_min_length_memory.get(&memory_key)
        {
            secondary_root_node.node_data.a_star_lower_bound += *a_star_lower_bound;
            context.open_list.clear();
            context.open_list.push(secondary_root_node);
        } else {
            context.open_list.clear();
            context.closed_list.clear();

            let mut open_list = mem::replace(&mut context.open_list, BinaryHeap::new_min());
            let mut closed_list = mem::take(&mut context.closed_list);
            let initial_cost = secondary_root_node.cost();
            open_list.push(secondary_root_node);

            let (target_node_identifier, _) = a_star_align_loop(
                reference,
                query,
                context,
                &mut closed_list,
                &mut open_list,
                |node, _, _, context| {
                    let Identifier::Secondary { length, .. } = node.node_data.identifier else {
                        unreachable!("A non-secondary node was closed before a target was closed.")
                    };
                    debug_assert!(length <= context.config.min_length);

                    length == context.config.min_length
                },
            );

            let target_cost = closed_list.get(&target_node_identifier).unwrap().cost();
            let lower_bound = target_cost - initial_cost;

            closed_nodes_output.extend(closed_list.drain().map(|(identifier, mut node)| {
                node.node_data.a_star_lower_bound += target_cost - node.node_data.cost;
                (identifier, node)
            }));

            context
                .template_switch_min_length_memory
                .insert(memory_key, lower_bound);

            context.open_list = open_list;
            context.closed_list = closed_list;
        }

        context.open_list.drain()
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
