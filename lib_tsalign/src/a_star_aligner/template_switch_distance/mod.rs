use compact_genome::interface::alphabet::AlphabetCharacter;
use identifier::{GapType, TemplateSwitchPrimary, TemplateSwitchSecondary};
use num_traits::SaturatingSub;
use strategies::{
    node_ord::NodeOrdStrategy, template_switch_min_length::TemplateSwitchMinLengthStrategy,
    AlignmentStrategies, AlignmentStrategy, AlignmentStrategySelector,
};

use crate::costs::cost::Cost;

use super::AlignmentGraphNode;

mod alignment_type;
mod context;
pub mod display;
mod identifier;
pub mod strategies;

pub use alignment_type::AlignmentType;
pub use context::Context;
pub use identifier::Identifier;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node<Strategies: AlignmentStrategySelector> {
    node_data: NodeData,
    strategies: AlignmentStrategies<Strategies>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct NodeData {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
    a_star_lower_bound: Cost,
}

impl<Strategies: AlignmentStrategySelector> AlignmentGraphNode<Strategies::Alphabet>
    for Node<Strategies>
{
    type Identifier = Identifier;

    type Context = Context<Strategies>;

    type AlignmentType = AlignmentType;

    fn create_root(context: &Self::Context) -> Self {
        Self {
            node_data: NodeData {
                identifier: Identifier::new_primary(0, 0, 0, GapType::None),
                predecessor: None,
                cost: Cost::ZERO,
                a_star_lower_bound: Cost::ZERO,
            },
            strategies: AlignmentStrategy::create_root(context),
        }
    }

    fn generate_successors<
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        context: &mut Self::Context,
        output: &mut impl Extend<Self>,
    ) {
        let config = &context.config;

        match self.node_data.identifier {
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
                debug_assert!(reference_index != usize::MAX, "{self:?}");
                debug_assert!(query_index != usize::MAX, "{self:?}");
                debug_assert!(reference_index < isize::MAX as usize, "{self:?}");
                debug_assert!(query_index < isize::MAX as usize, "{self:?}");

                if reference_index < reference.len() && query_index < query.len() {
                    // Diagonal characters
                    let r = reference[reference_index].clone();
                    let q = query[query_index].clone();

                    if flank_index == 0 {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .match_or_substitution_cost(r.clone(), q.clone()),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < config.left_flank_length {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .match_or_substitution_cost(r, q),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .match_or_substitution_cost(r, q),
                                context,
                            ),
                        );
                    } else if flank_index == config.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            config.offset_costs.evaluate(&0),
                            context,
                        ));
                    }
                }

                if reference_index < reference.len() {
                    // Deleted character
                    let r = reference[reference_index].clone();

                    if flank_index == 0 {
                        output.extend(
                            self.generate_primary_deletion_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .gap_costs(r.clone(), gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < config.left_flank_length {
                        output.extend(
                            self.generate_primary_deletion_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_deletion_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    } else if flank_index == config.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            config.offset_costs.evaluate(&0),
                            context,
                        ));
                    }
                }

                if query_index < query.len() {
                    // Inserted character
                    let q = query[query_index].clone();

                    if flank_index == 0 {
                        output.extend(
                            self.generate_primary_insertion_successor(
                                0,
                                config
                                    .primary_edit_costs
                                    .gap_costs(q.clone(), gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < config.left_flank_length {
                        output.extend(
                            self.generate_primary_insertion_successor(
                                flank_index + 1,
                                config
                                    .left_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_insertion_successor(
                                flank_index + 1,
                                config
                                    .right_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    } else if flank_index == config.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            config.offset_costs.evaluate(&0),
                            context,
                        ));
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
                let (secondary_entrance_index, secondary_length) = match template_switch_secondary {
                    TemplateSwitchSecondary::Reference => {
                        (entrance_reference_index, reference.len())
                    }
                    TemplateSwitchSecondary::Query => (entrance_query_index, query.len()),
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

                        output.extend(self.generate_template_switch_entrance_successor(
                            cost_increment,
                            template_switch_first_offset + 1,
                            context,
                        ))
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

                        output.extend(self.generate_template_switch_entrance_successor(
                            cost_increment,
                            template_switch_first_offset - 1,
                            context,
                        ))
                    }
                }

                output.extend(self.generate_secondary_root_node(reference, query, context));
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
                    TemplateSwitchPrimary::Reference => reference,
                    TemplateSwitchPrimary::Query => query,
                };
                let secondary_sequence = match template_switch_secondary {
                    TemplateSwitchSecondary::Reference => reference,
                    TemplateSwitchSecondary::Query => query,
                };

                if primary_index < primary_sequence.len() && secondary_index > 0 {
                    // Diagonal characters
                    let p = primary_sequence[primary_index].clone();
                    let s = secondary_sequence[secondary_index - 1].complement();

                    output.extend(self.generate_secondary_diagonal_successor(
                        config.secondary_edit_costs.match_or_substitution_cost(p, s),
                        context,
                    ));
                }

                if secondary_index > 0 {
                    // Deleted character
                    let s = secondary_sequence[secondary_index - 1].complement();

                    output.extend(
                        self.generate_secondary_deletion_successor(
                            config
                                .secondary_edit_costs
                                .gap_costs(s, gap_type != GapType::Deletion),
                            context,
                        ),
                    );
                }

                if primary_index < primary_sequence.len() {
                    // Inserted character
                    let p = primary_sequence[primary_index].clone();

                    output.extend(
                        self.generate_secondary_insertion_successor(
                            config
                                .secondary_edit_costs
                                .gap_costs(p, gap_type != GapType::Insertion),
                            context,
                        ),
                    );
                }

                let length_cost = config.length_costs.evaluate(&length);
                if length_cost != Cost::MAX {
                    let length_difference_cost = config.length_difference_costs.evaluate(&0);
                    assert_ne!(length_difference_cost, Cost::MAX);
                    let cost_increment = length_cost + length_difference_cost;

                    output.extend(
                        self.generate_initial_template_switch_exit_successor(
                            cost_increment,
                            context,
                        ),
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
                    TemplateSwitchPrimary::Reference => query.len(),
                    TemplateSwitchPrimary::Query => reference.len(),
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

                        output.extend(self.generate_template_switch_exit_successor(
                            cost_increment,
                            length_difference + 1,
                            context,
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

                        output.extend(self.generate_template_switch_exit_successor(
                            cost_increment,
                            length_difference - 1,
                            context,
                        ))
                    }
                }

                output.extend(self.generate_primary_reentry_successor(context));
            }
        }
    }

    fn identifier(&self) -> &Self::Identifier {
        &self.node_data.identifier
    }

    fn cost(&self) -> Cost {
        self.node_data.cost
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.node_data.predecessor.as_ref()
    }

    fn predecessor_alignment_type<
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        _context: &Self::Context,
    ) -> Self::AlignmentType {
        if self.node_data.predecessor.is_some() {
            match self.node_data.identifier {
                Identifier::Primary {
                    reference_index,
                    query_index,
                    gap_type,
                    ..
                } => match gap_type {
                    GapType::Insertion => AlignmentType::Insertion,
                    GapType::Deletion => AlignmentType::Deletion,
                    GapType::None => {
                        if reference[reference_index - 1] == query[query_index - 1] {
                            AlignmentType::Match
                        } else {
                            AlignmentType::Substitution
                        }
                    }
                },

                Identifier::PrimaryReentry { .. } => AlignmentType::PrimaryReentry,

                Identifier::TemplateSwitchEntrance {
                    template_switch_primary,
                    template_switch_secondary,
                    template_switch_first_offset,
                    ..
                } => AlignmentType::TemplateSwitchEntrance {
                    primary: template_switch_primary,
                    secondary: template_switch_secondary,
                    first_offset: template_switch_first_offset,
                },

                Identifier::Secondary {
                    template_switch_primary,
                    template_switch_secondary,
                    length,
                    primary_index,
                    secondary_index,
                    gap_type,
                    ..
                } => match gap_type {
                    GapType::Insertion => AlignmentType::Insertion,
                    GapType::Deletion => AlignmentType::Deletion,
                    GapType::None => {
                        if length == 0 {
                            AlignmentType::SecondaryRoot
                        } else {
                            let primary_character = match template_switch_primary {
                                TemplateSwitchPrimary::Reference => &reference[primary_index - 1],
                                TemplateSwitchPrimary::Query => &query[primary_index - 1],
                            };
                            let secondary_character = match template_switch_secondary {
                                TemplateSwitchSecondary::Reference => &reference[secondary_index],
                                TemplateSwitchSecondary::Query => &query[secondary_index],
                            };

                            if primary_character == &secondary_character.complement() {
                                AlignmentType::Match
                            } else {
                                AlignmentType::Substitution
                            }
                        }
                    }
                },

                Identifier::TemplateSwitchExit {
                    length_difference, ..
                } => AlignmentType::TemplateSwitchExit { length_difference },
            }
        } else {
            AlignmentType::Root
        }
    }

    fn is_target<
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        _context: &Self::Context,
    ) -> bool {
        match self.node_data.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                ..
            } => reference_index == reference.len() && query_index == query.len(),
            _ => false,
        }
    }
}

impl<Strategies: AlignmentStrategySelector> Node<Strategies> {
    fn generate_primary_diagonal_successor(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { .. }
        | Identifier::PrimaryReentry { .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_diagonal_successor(successor_flank_index),
            cost_increment,
            context,
        ))
    }

    fn generate_primary_deletion_successor(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { .. }
        | Identifier::PrimaryReentry { .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_deletion_successor(successor_flank_index),
            cost_increment,
            context,
        ))
    }

    fn generate_primary_insertion_successor(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ (Identifier::Primary { .. }
        | Identifier::PrimaryReentry { .. }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_insertion_successor(successor_flank_index),
            cost_increment,
            context,
        ))
    }

    fn generate_initial_template_switch_entrance_successors<'result>(
        &'result self,
        cost_increment: Cost,
        context: &'result <Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> impl 'result + Iterator<Item = Self> {
        if !matches!(
            self.node_data.identifier,
            Identifier::Primary { .. } | Identifier::PrimaryReentry { .. }
        ) {
            unreachable!("This method is only called on primary nodes.")
        }

        self.node_data
            .identifier
            .generate_initial_template_switch_entrance_successors()
            .map(move |identifier| self.generate_successor(identifier, cost_increment, context))
    }

    fn generate_template_switch_entrance_successor(
        &self,
        cost_increment: Cost,
        successor_template_switch_first_offset: isize,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_first_offset: successor_template_switch_first_offset,
            },
            cost_increment,
            context,
        ))
    }

    fn generate_secondary_root_node<
        'result,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &'result self,
        reference: &'result SubsequenceType,
        query: &'result SubsequenceType,
        context: &'result mut <Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> impl 'result + IntoIterator<Item = Self> {
        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_first_offset,
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        let primary_index = match template_switch_primary {
            TemplateSwitchPrimary::Reference => entrance_reference_index,
            TemplateSwitchPrimary::Query => entrance_query_index,
        };

        let secondary_index = (match template_switch_secondary {
            TemplateSwitchSecondary::Reference => entrance_reference_index,
            TemplateSwitchSecondary::Query => entrance_query_index,
        } as isize
            + template_switch_first_offset) as usize;

        let secondary_root_node = self.generate_successor(
            Identifier::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: 0,
                primary_index,
                secondary_index,
                gap_type: GapType::None,
            },
            0.into(),
            context,
        );

        self.strategies
            .template_switch_min_length_strategy
            .template_switch_min_length_lookahead(reference, query, secondary_root_node, context)
    }

    fn generate_secondary_diagonal_successor(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_diagonal_successor(),
            cost_increment,
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_deletion_successor(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_deletion_successor(),
            cost_increment,
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_insertion_successor(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_insertion_successor(),
            cost_increment,
            context,
        ))
    }

    fn generate_initial_template_switch_exit_successor(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::Secondary {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            primary_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                length_difference: 0,
            },
            cost_increment,
            context,
        ))
    }

    fn generate_template_switch_exit_successor(
        &self,
        cost_increment: Cost,
        successor_length_difference: isize,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            primary_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                length_difference: successor_length_difference,
            },
            cost_increment,
            context,
        ))
    }

    fn generate_primary_reentry_successor(
        &self,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        let Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            primary_index,
            length_difference,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        let (reference_index, query_index) = match template_switch_primary {
            TemplateSwitchPrimary::Reference => {
                let primary_length = primary_index - entrance_reference_index;
                let anti_primary_length = primary_length as isize + length_difference;
                let query_index = entrance_query_index as isize + anti_primary_length;

                if query_index < 0 {
                    return None;
                }

                (primary_index, query_index as usize)
            }
            TemplateSwitchPrimary::Query => {
                let primary_length = primary_index - entrance_query_index;
                let anti_primary_length = primary_length as isize + length_difference;
                let reference_index = entrance_reference_index as isize + anti_primary_length;

                if reference_index < 0 {
                    return None;
                }

                (reference_index as usize, primary_index)
            }
        };

        debug_assert!(reference_index != usize::MAX, "{self:?}");
        debug_assert!(query_index != usize::MAX, "{self:?}");
        debug_assert!(reference_index < isize::MAX as usize, "{self:?}");
        debug_assert!(query_index < isize::MAX as usize, "{self:?}");

        Some(self.generate_successor(
            Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: -context.config.right_flank_length,
            },
            0.into(),
            context,
        ))
    }

    fn generate_successor(
        &self,
        identifier: Identifier,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Self {
        Self {
            node_data: self
                .node_data
                .generate_successor(identifier, cost_increment),
            strategies: self.strategies.generate_successor(context),
        }
    }
}

impl NodeData {
    fn generate_successor(&self, identifier: Identifier, cost_increment: Cost) -> Self {
        let cost = self.cost + cost_increment;
        let a_star_lower_bound = self.a_star_lower_bound.saturating_sub(&cost_increment);
        Self {
            identifier,
            predecessor: Some(self.identifier),
            cost,
            a_star_lower_bound,
        }
    }
}

impl<Strategies: AlignmentStrategySelector> Ord for Node<Strategies> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.strategies
            .node_ord_strategy
            .cmp(&self.node_data, &other.node_data)
    }
}

impl<Strategies: AlignmentStrategySelector> PartialOrd for Node<Strategies> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
