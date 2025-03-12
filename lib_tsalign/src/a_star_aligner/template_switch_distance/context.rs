use std::fmt::Display;

use compact_genome::interface::alphabet::AlphabetCharacter;
use compact_genome::interface::sequence::GenomeSequence;
use extend_map::ExtendMap;
use generic_a_star::reset::Reset;
use generic_a_star::{AStarBuffers, AStarContext};
use num_traits::{Bounded, Zero};

use crate::a_star_aligner::template_switch_distance::Node;
use crate::a_star_aligner::AlignmentContext;
use crate::config::TemplateSwitchConfig;

use super::identifier::{GapType, TemplateSwitchPrimary, TemplateSwitchSecondary};
use super::strategies::chaining::ChainingStrategy;
use super::strategies::primary_match::PrimaryMatchStrategy;
use super::strategies::secondary_deletion::SecondaryDeletionStrategy;
use super::strategies::shortcut::ShortcutStrategy;
use super::strategies::template_switch_count::TemplateSwitchCountStrategy;
use super::strategies::template_switch_min_length::TemplateSwitchMinLengthStrategy;
use super::strategies::{AlignmentStrategiesNodeMemory, AlignmentStrategySelector};
use super::{AlignmentType, Identifier, NodeData};

pub struct Context<
    'reference,
    'query,
    SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    Strategies: AlignmentStrategySelector,
> {
    pub reference: &'reference SubsequenceType,
    pub query: &'query SubsequenceType,
    pub reference_name: String,
    pub query_name: String,

    pub config: TemplateSwitchConfig<Strategies::Alphabet, Strategies::Cost>,

    #[allow(clippy::type_complexity)]
    pub a_star_buffers: AStarBuffers<
        Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        Node<Strategies>,
    >,
    pub memory: Memory<Strategies>,

    cost_limit: Option<Strategies::Cost>,
    memory_limit: Option<usize>,
}

pub struct Memory<Strategies: AlignmentStrategySelector> {
    pub template_switch_min_length: <<Strategies as AlignmentStrategySelector>::TemplateSwitchMinLength as TemplateSwitchMinLengthStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::Memory,
    pub chaining: <<Strategies as AlignmentStrategySelector>::Chaining as ChainingStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::Memory,
    pub template_switch_count:  <<Strategies as AlignmentStrategySelector>::TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
    pub shortcut: <<Strategies as AlignmentStrategySelector>::Shortcut as ShortcutStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::Memory,
    pub primary_match: <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
    <Strategies as AlignmentStrategySelector>::Cost,
>>::Memory,
}

impl<
        'reference,
        'query,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > Context<'reference, 'query, SubsequenceType, Strategies>
{
    #[expect(clippy::too_many_arguments)]
    pub fn new(
        reference: &'reference SubsequenceType,
        query: &'query SubsequenceType,
        reference_name: &str,
        query_name: &str,
        config: TemplateSwitchConfig<Strategies::Alphabet, Strategies::Cost>,
        memory: Memory<Strategies>,
        cost_limit: Option<Strategies::Cost>,
        memory_limit: Option<usize>,
    ) -> Self {
        Self {
            reference,
            query,
            reference_name: reference_name.to_owned(),
            query_name: query_name.to_owned(),
            config,
            a_star_buffers: Default::default(),
            memory,
            cost_limit,
            memory_limit,
        }
    }
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > AStarContext for Context<'_, '_, SubsequenceType, Strategies>
{
    type Node = Node<Strategies>;

    fn create_root(&self) -> Self::Node {
        Self::Node {
            node_data: NodeData {
                identifier: Identifier::new_primary(0, 0, 0, GapType::None, <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::create_root_identifier_primary_extra_data(self)),
                predecessor: None,
                predecessor_edge_type: AlignmentType::Root,
                cost: Strategies::Cost::zero(),
                a_star_lower_bound: Strategies::Cost::zero(),
            },
            strategies: AlignmentStrategiesNodeMemory::create_root(self),
        }
    }

    fn generate_successors(
        &mut self,
        node: &Self::Node,
        opened_nodes_output: &mut impl Extend<Self::Node>,
    ) {
        let config = &self.config;
        let mut opened_nodes_output =
            ExtendMap::new(opened_nodes_output, generate_output_mapper_function(self));

        match node.node_data.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
                ..
            }
            | Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type,
                flank_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX, "{node:?}");
                debug_assert!(query_index != usize::MAX, "{node:?}");
                debug_assert!(reference_index < isize::MAX as usize, "{node:?}");
                debug_assert!(query_index < isize::MAX as usize, "{node:?}");

                let can_start_another_template_switch = node
                    .strategies
                    .template_switch_count
                    .can_start_another_template_switch(self);

                if reference_index < self.reference.len() && query_index < self.query.len() {
                    // Diagonal characters
                    let r = self.reference[reference_index].clone();
                    let q = self.query[query_index].clone();
                    let is_match = r == q;

                    if flank_index == 0 {
                        let can_do_primary_non_flank_match = <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::can_do_primary_non_flank_match(node.node_data.identifier, self);

                        let (is_match, cost_increment) =
                            if is_match && can_do_primary_non_flank_match {
                                (
                                    true,
                                    config.primary_edit_costs.match_cost(r.clone(), q.clone()),
                                )
                            } else if is_match && !can_do_primary_non_flank_match {
                                (
                                    false,
                                    node.strategies.primary_match.fake_substitution_cost(self),
                                )
                            } else {
                                debug_assert!(!is_match);
                                (
                                    false,
                                    config
                                        .primary_edit_costs
                                        .substitution_cost(r.clone(), q.clone()),
                                )
                            };

                        if cost_increment != Strategies::Cost::max_value() {
                            opened_nodes_output.extend(node.generate_primary_diagonal_successor(
                                0,
                                cost_increment,
                                is_match,
                                self,
                            ));
                        }

                        if is_match && <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::always_generate_substitution() {
                            let cost_increment = node.strategies.primary_match.fake_substitution_cost(self);

                            if cost_increment != Strategies::Cost::max_value() {
                                opened_nodes_output.extend(node.generate_primary_diagonal_successor(
                                    0,
                                    cost_increment,
                                    false,
                                    self,
                                ));
                            }
                        }
                    }

                    if (flank_index < config.left_flank_length && can_start_another_template_switch)
                        || flank_index < 0
                    {
                        let can_do_primary_flank_match = <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<<Strategies as AlignmentStrategySelector>::Cost>>::can_do_primary_flank_match(node.node_data.identifier, self);

                        let edit_costs = if flank_index < 0 {
                            &config.right_flank_edit_costs
                        } else {
                            &config.left_flank_edit_costs
                        };

                        let (is_match, cost_increment) = if is_match && can_do_primary_flank_match {
                            (true, edit_costs.match_cost(r.clone(), q.clone()))
                        } else if is_match && !can_do_primary_flank_match {
                            (
                                false,
                                node.strategies.primary_match.fake_substitution_cost(self),
                            )
                        } else {
                            debug_assert!(!is_match);
                            (false, edit_costs.substitution_cost(r.clone(), q.clone()))
                        };

                        if cost_increment != Strategies::Cost::max_value() {
                            opened_nodes_output.extend(node.generate_primary_diagonal_successor(
                                flank_index + 1,
                                cost_increment,
                                is_match,
                                self,
                            ));
                        }
                    }
                }

                if reference_index < self.reference.len() {
                    // Deleted character
                    let r = self.reference[reference_index].clone();

                    if flank_index == 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_deletion_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .gap_costs(r.clone(), gap_type != GapType::Deletion),
                                self,
                            ),
                        );
                    }

                    if flank_index >= 0
                        && flank_index < config.left_flank_length
                        && can_start_another_template_switch
                    {
                        opened_nodes_output.extend(
                            node.generate_primary_deletion_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                self,
                            ),
                        );
                    } else if flank_index < 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_deletion_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                self,
                            ),
                        );
                    }
                }

                if query_index < self.query.len() {
                    // Inserted character
                    let q = self.query[query_index].clone();

                    if flank_index == 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_insertion_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .gap_costs(q.clone(), gap_type != GapType::Insertion),
                                self,
                            ),
                        );
                    }

                    if flank_index >= 0
                        && flank_index < config.left_flank_length
                        && can_start_another_template_switch
                    {
                        opened_nodes_output.extend(
                            node.generate_primary_insertion_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                self,
                            ),
                        );
                    } else if flank_index < 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_insertion_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                self,
                            ),
                        );
                    }
                }

                // Template switches are always allowed, as long as we have a left flank.
                if flank_index == config.left_flank_length && can_start_another_template_switch {
                    let offset_costs = config.offset_costs.evaluate(&0);

                    if offset_costs != Strategies::Cost::max_value()
                        && config.base_cost != Strategies::Cost::max_value()
                    {
                        opened_nodes_output.extend(
                            node.generate_initial_template_switch_entrance_successors(
                                config.offset_costs.evaluate(&0) + config.base_cost,
                                self,
                            ),
                        );
                    }
                }
            }

            Identifier::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_secondary,
                template_switch_first_offset,
                ..
            } => {
                debug_assert!(node
                    .strategies
                    .template_switch_count
                    .can_start_another_template_switch(self));

                let (secondary_entrance_index, secondary_length) = match template_switch_secondary {
                    TemplateSwitchSecondary::Reference => {
                        (entrance_reference_index, self.reference.len())
                    }
                    TemplateSwitchSecondary::Query => (entrance_query_index, self.query.len()),
                };
                let secondary_index =
                    secondary_entrance_index as isize + template_switch_first_offset;

                if template_switch_first_offset >= 0 && secondary_index < secondary_length as isize
                {
                    let new_cost = config
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset + 1));

                    if new_cost != Strategies::Cost::max_value() {
                        let old_cost = config.offset_costs.evaluate(&template_switch_first_offset);
                        assert!(new_cost >= old_cost);
                        let cost_increment = new_cost - old_cost;

                        opened_nodes_output.extend(
                            node.generate_template_switch_entrance_successor(
                                cost_increment,
                                template_switch_first_offset + 1,
                                self,
                            ),
                        )
                    }
                }

                if template_switch_first_offset <= 0
                    && secondary_index > self.config.min_length as isize
                {
                    let new_cost = config
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset - 1));

                    if new_cost != Strategies::Cost::max_value() {
                        let old_cost = config.offset_costs.evaluate(&template_switch_first_offset);
                        assert!(new_cost >= old_cost);
                        let cost_increment = new_cost - old_cost;

                        opened_nodes_output.extend(
                            node.generate_template_switch_entrance_successor(
                                cost_increment,
                                template_switch_first_offset - 1,
                                self,
                            ),
                        )
                    }
                }

                if secondary_index >= self.config.min_length as isize
                    && secondary_index <= secondary_length as isize
                {
                    // Temporarily unpack opened_nodes_output because it borrows self,
                    // but generating the secondary root node wants to borrow self as mutable.
                    let opened_nodes_direct_output = opened_nodes_output.into_inner();
                    let secondary_root_node: Vec<_> = node
                        .generate_secondary_root_node(self)
                        .into_iter()
                        .collect();
                    opened_nodes_output = ExtendMap::new(
                        opened_nodes_direct_output,
                        generate_output_mapper_function(self),
                    );
                    opened_nodes_output.extend(secondary_root_node);
                }
            }

            Identifier::Secondary {
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                gap_type,
                ..
            } => {
                // TODO Some of the nodes generated here are unable to reach the target:
                // TODO * nodes who get closer than `right_flank_length` to the end of the primary sequence
                // TODO * nodes who get closer than `right_flank_length` to the end of the not-primary sequence,
                // TODO   assuming a `length_difference` of zero

                let primary_sequence = match template_switch_primary {
                    TemplateSwitchPrimary::Reference => self.reference,
                    TemplateSwitchPrimary::Query => self.query,
                };
                let secondary_sequence = match template_switch_secondary {
                    TemplateSwitchSecondary::Reference => self.reference,
                    TemplateSwitchSecondary::Query => self.query,
                };

                // Only generate secondary successors if they can ever exit the template switch based on their length.
                let min_length_cost = config.length_costs.min(length..).unwrap();
                if min_length_cost != Strategies::Cost::max_value() {
                    if primary_index < primary_sequence.len() && secondary_index > 0 {
                        // Diagonal characters
                        let p = primary_sequence[primary_index].clone();
                        let s = secondary_sequence[secondary_index - 1].complement();

                        opened_nodes_output.extend(
                            node.generate_secondary_diagonal_successor(
                                config
                                    .secondary_edit_costs
                                    .match_or_substitution_cost(p.clone(), s.clone()),
                                p == s,
                                self,
                            ),
                        );
                    }

                    if secondary_index > 0
                        && Strategies::SecondaryDeletion::allow_secondary_deletions()
                    {
                        if secondary_index > secondary_sequence.len() {
                            panic!("Secondary index out of bounds for node {node}");
                        }

                        // Deleted character
                        let s = secondary_sequence[secondary_index - 1].complement();

                        opened_nodes_output.extend(
                            node.generate_secondary_deletion_successor(
                                config
                                    .secondary_edit_costs
                                    .gap_costs(s, gap_type != GapType::Deletion),
                                self,
                            ),
                        );
                    }

                    if primary_index < primary_sequence.len() {
                        // Inserted character
                        let p = primary_sequence[primary_index].clone();

                        opened_nodes_output.extend(
                            node.generate_secondary_insertion_successor(
                                config
                                    .secondary_edit_costs
                                    .gap_costs(p, gap_type != GapType::Insertion),
                                self,
                            ),
                        );
                    }
                }

                let length_cost = config.length_costs.evaluate(&length);
                if length_cost != Strategies::Cost::max_value() {
                    let length_difference_cost = config.length_difference_costs.evaluate(&0);
                    assert_ne!(length_difference_cost, Strategies::Cost::max_value());
                    let cost_increment = length_cost + length_difference_cost;

                    opened_nodes_output.extend(
                        node.generate_initial_template_switch_exit_successor(cost_increment, self),
                    )
                }
            }

            Identifier::TemplateSwitchExit {
                template_switch_primary,
                primary_index,
                length_difference,
                ..
            } => {
                let anti_primary_length = match template_switch_primary {
                    TemplateSwitchPrimary::Reference => self.query.len(),
                    TemplateSwitchPrimary::Query => self.reference.len(),
                };

                if length_difference >= 0
                    && primary_index as isize + length_difference < anti_primary_length as isize
                {
                    let new_cost = config
                        .length_difference_costs
                        .evaluate(&(&length_difference + 1));

                    if new_cost != Strategies::Cost::max_value() {
                        let old_cost = config.length_difference_costs.evaluate(&length_difference);
                        assert!(new_cost >= old_cost);
                        let cost_increment = new_cost - old_cost;

                        opened_nodes_output.extend(node.generate_template_switch_exit_successor(
                            cost_increment,
                            length_difference + 1,
                            self,
                        ))
                    }
                }

                if length_difference <= 0 && primary_index as isize + length_difference > 0 {
                    let new_cost = config
                        .length_difference_costs
                        .evaluate(&(&length_difference - 1));

                    if new_cost != Strategies::Cost::max_value() {
                        let old_cost = config.length_difference_costs.evaluate(&length_difference);
                        assert!(new_cost >= old_cost);
                        let cost_increment = new_cost - old_cost;

                        opened_nodes_output.extend(node.generate_template_switch_exit_successor(
                            cost_increment,
                            length_difference - 1,
                            self,
                        ))
                    }
                }

                opened_nodes_output.extend(node.generate_primary_reentry_successor(self).map(
                    |mut node| {
                        node.strategies.template_switch_count.increment_count();
                        node
                    },
                ));
            }
        }

        // Add additional successors through strategies.
        <<Strategies as AlignmentStrategySelector>::Shortcut as ShortcutStrategy<
            <Strategies as AlignmentStrategySelector>::Cost,
        >>::generate_successors(node, self, &mut opened_nodes_output);
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        match node.node_data.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                ..
            } => reference_index == self.reference.len() && query_index == self.query.len(),
            _ => false,
        }
    }

    fn cost_limit(&self) -> Option<Strategies::Cost> {
        self.cost_limit
    }

    fn memory_limit(&self) -> Option<usize> {
        self.memory_limit
    }
}

fn generate_output_mapper_function<
    'context,
    'reference,
    'query,
    SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    Strategies: AlignmentStrategySelector,
>(
    context: &'context Context<'reference, 'query, SubsequenceType, Strategies>,
) -> impl use<'context, 'reference, 'query, SubsequenceType, Strategies>
       + Fn(
    <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node,
) -> <Context<'reference, 'query, SubsequenceType, Strategies> as AStarContext>::Node {
    move |node| {
        <Strategies as AlignmentStrategySelector>::Chaining::apply_lower_bound(node, context)
    }
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > Reset for Context<'_, '_, SubsequenceType, Strategies>
{
    fn reset(&mut self) {
        self.memory.reset();
    }
}

impl<Strategies: AlignmentStrategySelector> Reset for Memory<Strategies> {
    fn reset(&mut self) {
        self.template_switch_min_length.reset();
    }
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > AlignmentContext for Context<'_, '_, SubsequenceType, Strategies>
{
    type AlphabetType = Strategies::Alphabet;

    type AlignmentType = AlignmentType;

    type SubsequenceType = SubsequenceType;

    fn reference(&self) -> &Self::SubsequenceType {
        self.reference
    }

    fn query(&self) -> &Self::SubsequenceType {
        self.query
    }

    fn reference_name(&self) -> &str {
        &self.reference_name
    }

    fn query_name(&self) -> &str {
        &self.query_name
    }
}

impl<Strategies: AlignmentStrategySelector> Display for AlignmentStrategiesNodeMemory<Strategies> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.primary_match)
    }
}
