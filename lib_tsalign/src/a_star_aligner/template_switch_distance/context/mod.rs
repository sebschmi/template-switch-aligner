use compact_genome::interface::alphabet::AlphabetCharacter;
use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::cost::Cost;
use generic_a_star::reset::Reset;
use generic_a_star::{AStarBuffers, AStarContext};
use log::debug;

use crate::a_star_aligner::template_switch_distance::Node;
use crate::a_star_aligner::AlignmentContext;
use crate::config::TemplateSwitchConfig;

use super::identifier::{GapType, TemplateSwitchPrimary, TemplateSwitchSecondary};
use super::strategies::chaining::ChainingStrategy;
use super::strategies::secondary_deletion_strategy::SecondaryDeletionStrategy;
use super::strategies::shortcut::ShortcutStrategy;
use super::strategies::template_switch_count::TemplateSwitchCountStrategy;
use super::strategies::template_switch_min_length::TemplateSwitchMinLengthStrategy;
use super::strategies::{AlignmentStrategy, AlignmentStrategySelector};
use super::{AlignmentType, Identifier, NodeData};

pub struct Context<'reference, 'query, SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized, Strategies: AlignmentStrategySelector> {
    pub reference: &'reference SubsequenceType,
    pub query: &'query SubsequenceType,

    pub config: TemplateSwitchConfig<Strategies::Alphabet>,

    pub template_switch_min_length_memory: <<Strategies as AlignmentStrategySelector>::TemplateSwitchMinLength as TemplateSwitchMinLengthStrategy>::Memory,
    pub a_star_buffers: AStarBuffers<Identifier, Node<Strategies>>,
    pub chaining_memory: <<Strategies as AlignmentStrategySelector>::Chaining as ChainingStrategy>::Memory,
    pub template_switch_count_memory:  <<Strategies as AlignmentStrategySelector>::TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
}

impl<
        'reference,
        'query,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > Context<'reference, 'query, SubsequenceType, Strategies>
{
    pub fn new(
        reference: &'reference SubsequenceType,
        query: &'query SubsequenceType,
        config: TemplateSwitchConfig<Strategies::Alphabet>,
        template_switch_count_memory:<<Strategies as AlignmentStrategySelector>::TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
    ) -> Self {
        debug!("Creating/loading context...");
        let chaining_memory = <<Strategies as AlignmentStrategySelector>::Chaining as ChainingStrategy>::initialise_memory(&config);

        Self {
            reference,
            query,
            config,
            template_switch_min_length_memory: Default::default(),
            a_star_buffers: Default::default(),
            chaining_memory,
            template_switch_count_memory,
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
                identifier: Identifier::new_primary(0, 0, 0, GapType::None),
                predecessor: None,
                predecessor_edge_type: AlignmentType::Root,
                cost: Cost::ZERO,
                a_star_lower_bound: Cost::ZERO,
            },
            strategies: AlignmentStrategy::create_root(self),
        }
    }

    fn generate_successors(
        &mut self,
        node: &Self::Node,
        opened_nodes_output: &mut impl Extend<Self::Node>,
    ) {
        let config = &self.config;

        match node.node_data.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
            }
            | Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type,
                flank_index,
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

                    if flank_index == 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_diagonal_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .match_or_substitution_cost(r.clone(), q.clone()),
                                r == q,
                                self,
                            ),
                        );
                    }

                    if flank_index >= 0
                        && flank_index < config.left_flank_length
                        && can_start_another_template_switch
                    {
                        opened_nodes_output.extend(
                            node.generate_primary_diagonal_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .match_or_substitution_cost(r.clone(), q.clone()),
                                r == q,
                                self,
                            ),
                        );
                    } else if flank_index < 0 {
                        opened_nodes_output.extend(
                            node.generate_primary_diagonal_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .match_or_substitution_cost(r.clone(), q.clone()),
                                r == q,
                                self,
                            ),
                        );
                    } else if flank_index == config.left_flank_length
                        && can_start_another_template_switch
                    {
                        opened_nodes_output.extend(
                            node.generate_initial_template_switch_entrance_successors(
                                config.offset_costs.evaluate(&0) + config.base_cost,
                                self,
                            ),
                        );
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
                    } else if flank_index == config.left_flank_length
                        && can_start_another_template_switch
                    {
                        opened_nodes_output.extend(
                            node.generate_initial_template_switch_entrance_successors(
                                config.offset_costs.evaluate(&0) + config.base_cost,
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
                    } else if flank_index == config.left_flank_length
                        && can_start_another_template_switch
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

                if template_switch_first_offset >= 0
                    && secondary_entrance_index as isize + template_switch_first_offset
                        < secondary_length as isize
                {
                    let new_cost = config
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset + 1));

                    if new_cost != Cost::MAX {
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
                    && secondary_entrance_index as isize + template_switch_first_offset > 0
                {
                    let new_cost = config
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset - 1));

                    if new_cost != Cost::MAX {
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

                opened_nodes_output.extend(node.generate_secondary_root_node(self));
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
                if min_length_cost != Cost::MAX {
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
                if length_cost != Cost::MAX {
                    let length_difference_cost = config.length_difference_costs.evaluate(&0);
                    assert_ne!(length_difference_cost, Cost::MAX);
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

                    if new_cost != Cost::MAX {
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

                    if new_cost != Cost::MAX {
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
        <<Strategies as AlignmentStrategySelector>::Shortcut as ShortcutStrategy>::generate_successors(node, self, opened_nodes_output);
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
}

impl<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    > Reset for Context<'_, '_, SubsequenceType, Strategies>
{
    fn reset(&mut self) {
        self.template_switch_min_length_memory.reset();
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
}
