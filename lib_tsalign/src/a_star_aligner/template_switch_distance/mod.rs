use compact_genome::interface::alphabet::AlphabetCharacter;
use strategies::{
    node_ord::NodeOrdStrategy, AlignmentStrategies, AlignmentStrategy, AlignmentStrategySelector,
};

use crate::{config::TemplateSwitchConfig, costs::cost::Cost};

use super::{alignment_result::IAlignmentType, AlignmentGraphNode};

pub mod display;
pub mod strategies;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node<Strategies: AlignmentStrategySelector> {
    node_data: NodeData,
    strategies: AlignmentStrategies<Strategies>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NodeData {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Identifier {
    Primary {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
        /// Positive for left flank, negative for right flank.
        flank_index: isize,
    },
    TemplateSwitchEntrance {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        template_switch_first_offset: isize,
    },
    Secondary {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        length: usize,
        /// The index that does not jump.
        primary_index: usize,
        /// The index that jumps.
        secondary_index: usize,
        gap_type: GapType,
    },
    TemplateSwitchExit {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        /// The index that does not jump.
        primary_index: usize,
        length_difference: isize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GapType {
    Insertion,
    Deletion,
    None,
}

/// The primary sequence is the sequence for which the template switch does not jump.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TemplateSwitchPrimary {
    Reference,
    Query,
}

/// The secondary sequence is the sequence for which the template switch jumps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TemplateSwitchSecondary {
    Reference,
    Query,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    /// The query contains a base that is missing from the reference.
    Insertion,
    /// The query is missing a base present in the reference.
    Deletion,
    /// The query contains a different base than the reference.
    Substitution,
    /// The query contains the same base as the reference.
    Match,
    /// A template switch entrance.
    TemplateSwitchEntrance {
        primary: TemplateSwitchPrimary,
        secondary: TemplateSwitchSecondary,
        first_offset: isize,
    },
    /// A template switch exit.
    TemplateSwitchExit { length_difference: isize },
    /// This node is the root node, hence it was not generated via alignment.
    Root,
    /// The root node of a secondary graph.
    SecondaryRoot,
}

impl<Strategies: AlignmentStrategySelector> AlignmentGraphNode<Strategies::Alphabet>
    for Node<Strategies>
{
    type Identifier = Identifier;

    type Context = TemplateSwitchConfig<Strategies::Alphabet>;

    type AlignmentType = AlignmentType;

    fn create_root(context: &Self::Context) -> Self {
        Self {
            node_data: NodeData {
                identifier: Identifier::new_primary(0, 0, 0, GapType::None),
                predecessor: None,
                cost: Cost::ZERO,
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
        context: &Self::Context,
        output: &mut impl Extend<Self>,
    ) {
        match self.node_data.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
            } => {
                if reference_index < reference.len() && query_index < query.len() {
                    // Diagonal characters
                    let r = reference[reference_index].clone();
                    let q = query[query_index].clone();

                    if flank_index == 0 {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                0,
                                context
                                    .primary_edit_costs
                                    .match_or_substitution_cost(r.clone(), q.clone()),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                flank_index + 1,
                                context
                                    .left_flank_edit_costs
                                    .match_or_substitution_cost(r, q),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_diagonal_successor(
                                flank_index + 1,
                                context
                                    .right_flank_edit_costs
                                    .match_or_substitution_cost(r, q),
                                context,
                            ),
                        );
                    } else if flank_index == context.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            context.offset_costs.evaluate(&0),
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
                                context
                                    .primary_edit_costs
                                    .gap_costs(r.clone(), gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend(
                            self.generate_primary_deletion_successor(
                                flank_index + 1,
                                context
                                    .left_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_deletion_successor(
                                flank_index + 1,
                                context
                                    .right_flank_edit_costs
                                    .gap_costs(r, gap_type != GapType::Deletion),
                                context,
                            ),
                        );
                    } else if flank_index == context.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            context.offset_costs.evaluate(&0),
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
                                context
                                    .primary_edit_costs
                                    .gap_costs(q.clone(), gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend(
                            self.generate_primary_insertion_successor(
                                flank_index + 1,
                                context
                                    .left_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    } else if flank_index < 0 {
                        output.extend(
                            self.generate_primary_insertion_successor(
                                flank_index + 1,
                                context
                                    .right_flank_edit_costs
                                    .gap_costs(q, gap_type != GapType::Insertion),
                                context,
                            ),
                        );
                    } else if flank_index == context.left_flank_length {
                        output.extend(self.generate_initial_template_switch_entrance_successors(
                            context.offset_costs.evaluate(&0),
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
                    let new_cost = context
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset + 1));

                    if new_cost != Cost::MAX {
                        let old_cost = context.offset_costs.evaluate(&template_switch_first_offset);
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
                    let new_cost = context
                        .offset_costs
                        .evaluate(&(&template_switch_first_offset - 1));

                    if new_cost != Cost::MAX {
                        let old_cost = context.offset_costs.evaluate(&template_switch_first_offset);
                        assert!(new_cost >= old_cost);
                        let cost_increment = new_cost - old_cost;

                        output.extend(self.generate_template_switch_entrance_successor(
                            cost_increment,
                            template_switch_first_offset - 1,
                            context,
                        ))
                    }
                }

                output.extend(self.generate_secondary_root_node(context));
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
                    let s = secondary_sequence[secondary_index - 1].clone();

                    output.extend(
                        self.generate_secondary_diagonal_successor(
                            context
                                .secondary_edit_costs
                                .match_or_substitution_cost(p, s),
                            context,
                        ),
                    );
                }

                if secondary_index > 0 {
                    // Deleted character
                    let s = secondary_sequence[secondary_index - 1].clone();

                    output.extend(
                        self.generate_secondary_deletion_successor(
                            context
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
                            context
                                .secondary_edit_costs
                                .gap_costs(p, gap_type != GapType::Insertion),
                            context,
                        ),
                    );
                }

                let length_cost = context.length_costs.evaluate(&length);
                if length_cost != Cost::MAX {
                    let length_difference_cost = context.length_difference_costs.evaluate(&0);
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

            #[expect(unused)]
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                length_difference,
            } => {
                let anti_primary_length = match template_switch_primary {
                    TemplateSwitchPrimary::Reference => query.len(),
                    TemplateSwitchPrimary::Query => reference.len(),
                };

                if length_difference >= 0
                    && primary_index as isize + length_difference < anti_primary_length as isize
                {
                    let new_cost = context
                        .length_difference_costs
                        .evaluate(&(&length_difference + 1));

                    if new_cost != Cost::MAX {
                        let old_cost = context.length_difference_costs.evaluate(&length_difference);
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
                    let new_cost = context
                        .length_difference_costs
                        .evaluate(&(&length_difference - 1));

                    if new_cost != Cost::MAX {
                        let old_cost = context.length_difference_costs.evaluate(&length_difference);
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

        let predecessor_identifier @ Identifier::Primary { .. } = self.node_data.identifier else {
            unreachable!("This method is only called on primary nodes.")
        };

        let cost = self.node_data.cost + cost_increment;
        let node_data = NodeData {
            identifier: predecessor_identifier
                .generate_primary_diagonal_successor(successor_flank_index),
            predecessor: Some(predecessor_identifier),
            cost,
        };
        Some(Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        })
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

        let predecessor_identifier @ Identifier::Primary { .. } = self.node_data.identifier else {
            unreachable!("This method is only called on primary nodes.")
        };

        let cost = self.node_data.cost + cost_increment;
        let node_data = NodeData {
            identifier: predecessor_identifier
                .generate_primary_deletion_successor(successor_flank_index),
            predecessor: Some(predecessor_identifier),
            cost,
        };
        Some(Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        })
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

        let predecessor_identifier @ Identifier::Primary { .. } = self.node_data.identifier else {
            unreachable!("This method is only called on primary nodes.")
        };

        let cost = self.node_data.cost + cost_increment;
        let node_data = NodeData {
            identifier: predecessor_identifier
                .generate_primary_insertion_successor(successor_flank_index),
            predecessor: Some(predecessor_identifier),
            cost,
        };
        Some(Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        })
    }

    fn generate_initial_template_switch_entrance_successors(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> impl IntoIterator<Item = Self> {
        let predecessor_identifier @ Identifier::Primary {
            reference_index: entrance_reference_index,
            query_index: entrance_query_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        let template_switch_first_offset = 0;
        let predecessor = Some(predecessor_identifier);
        let cost = self.node_data.cost + cost_increment;

        [
            Self {
                node_data: NodeData {
                    identifier: Identifier::TemplateSwitchEntrance {
                        entrance_reference_index,
                        entrance_query_index,
                        template_switch_primary: TemplateSwitchPrimary::Reference,
                        template_switch_secondary: TemplateSwitchSecondary::Reference,
                        template_switch_first_offset,
                    },
                    predecessor,
                    cost,
                },
                strategies: self.strategies.generate_successor(context),
            },
            Self {
                node_data: NodeData {
                    identifier: Identifier::TemplateSwitchEntrance {
                        entrance_reference_index,
                        entrance_query_index,
                        template_switch_primary: TemplateSwitchPrimary::Reference,
                        template_switch_secondary: TemplateSwitchSecondary::Query,
                        template_switch_first_offset,
                    },
                    predecessor,
                    cost,
                },
                strategies: self.strategies.generate_successor(context),
            },
            Self {
                node_data: NodeData {
                    identifier: Identifier::TemplateSwitchEntrance {
                        entrance_reference_index,
                        entrance_query_index,
                        template_switch_primary: TemplateSwitchPrimary::Query,
                        template_switch_secondary: TemplateSwitchSecondary::Reference,
                        template_switch_first_offset,
                    },
                    predecessor,
                    cost,
                },
                strategies: self.strategies.generate_successor(context),
            },
            Self {
                node_data: NodeData {
                    identifier: Identifier::TemplateSwitchEntrance {
                        entrance_reference_index,
                        entrance_query_index,
                        template_switch_primary: TemplateSwitchPrimary::Query,
                        template_switch_secondary: TemplateSwitchSecondary::Query,
                        template_switch_first_offset,
                    },
                    predecessor,
                    cost,
                },
                strategies: self.strategies.generate_successor(context),
            },
        ]
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

        let predecessor_identifier @ Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        Some(Self {
            node_data: NodeData {
                identifier: Identifier::TemplateSwitchEntrance {
                    entrance_reference_index,
                    entrance_query_index,
                    template_switch_primary,
                    template_switch_secondary,
                    template_switch_first_offset: successor_template_switch_first_offset,
                },
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
    }

    fn generate_secondary_root_node(
        &self,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        let predecessor_identifier @ Identifier::TemplateSwitchEntrance {
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

        Some(Self {
            node_data: NodeData {
                identifier: Identifier::Secondary {
                    entrance_reference_index,
                    entrance_query_index,
                    template_switch_primary,
                    template_switch_secondary,
                    length: 0,
                    primary_index,
                    secondary_index,
                    gap_type: GapType::None,
                },
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost,
            },
            strategies: self.strategies.generate_successor(context),
        })
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

        Some(Self {
            node_data: NodeData {
                identifier: predecessor_identifier.generate_secondary_diagonal_successor(),
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
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

        Some(Self {
            node_data: NodeData {
                identifier: predecessor_identifier.generate_secondary_deletion_successor(),
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
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

        Some(Self {
            node_data: NodeData {
                identifier: predecessor_identifier.generate_secondary_insertion_successor(),
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
    }

    fn generate_initial_template_switch_exit_successor(
        &self,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary {
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

        Some(Self {
            node_data: NodeData {
                identifier: Identifier::TemplateSwitchExit {
                    entrance_reference_index,
                    entrance_query_index,
                    template_switch_primary,
                    template_switch_secondary,
                    primary_index,
                    length_difference: 0,
                },
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
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

        let predecessor_identifier @ Identifier::TemplateSwitchExit {
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

        Some(Self {
            node_data: NodeData {
                identifier: Identifier::TemplateSwitchExit {
                    entrance_reference_index,
                    entrance_query_index,
                    template_switch_primary,
                    template_switch_secondary,
                    primary_index,
                    length_difference: successor_length_difference,
                },
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost + cost_increment,
            },
            strategies: self.strategies.generate_successor(context),
        })
    }

    fn generate_primary_reentry_successor(
        &self,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Option<Self> {
        let predecessor_identifier @ Identifier::TemplateSwitchExit {
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
                let anti_primary_length = (primary_length as isize + length_difference) as usize;
                (primary_index, entrance_query_index + anti_primary_length)
            }
            TemplateSwitchPrimary::Query => {
                let primary_length = primary_index - entrance_query_index;
                let anti_primary_length = (primary_length as isize + length_difference) as usize;
                (
                    entrance_reference_index + anti_primary_length,
                    primary_index,
                )
            }
        };

        Some(Self {
            node_data: NodeData {
                identifier: Identifier::Primary {
                    reference_index,
                    query_index,
                    gap_type: GapType::None,
                    flank_index: -context.right_flank_length,
                },
                predecessor: Some(predecessor_identifier),
                cost: self.node_data.cost,
            },
            strategies: self.strategies.generate_successor(context),
        })
    }
}

impl Identifier {
    const fn new_primary(
        reference_index: usize,
        query_index: usize,
        flank_index: isize,
        gap_type: GapType,
    ) -> Self {
        Self::Primary {
            reference_index,
            query_index,
            flank_index,
            gap_type,
        }
    }

    fn generate_primary_diagonal_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index: query_index + 1,
                flank_index,
                gap_type: GapType::None,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    fn generate_primary_deletion_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index,
                flank_index,
                gap_type: GapType::Deletion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    fn generate_primary_insertion_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index,
                query_index: query_index + 1,
                flank_index,
                gap_type: GapType::Insertion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    fn generate_secondary_diagonal_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: length + 1,
                primary_index: primary_index + 1,
                secondary_index: secondary_index - 1,
                gap_type: GapType::None,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_deletion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index: secondary_index - 1,
                gap_type: GapType::Deletion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The secondary misses a base present in the primary.
    fn generate_secondary_insertion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: length + 1,
                primary_index: primary_index + 1,
                secondary_index,
                gap_type: GapType::Insertion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// Returns the anti-diagonal for variants where it exists, or [`usize::MAX`](core::primitive::usize::MAX) otherwise.
    const fn anti_diagonal(self) -> usize {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => reference_index + query_index,
            _ => usize::MAX,
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

impl IAlignmentType for AlignmentType {
    fn is_repeatable(&self) -> bool {
        match self {
            AlignmentType::Insertion
            | AlignmentType::Deletion
            | AlignmentType::Substitution
            | AlignmentType::Match
            | AlignmentType::Root
            | AlignmentType::SecondaryRoot => true,
            AlignmentType::TemplateSwitchEntrance { .. }
            | AlignmentType::TemplateSwitchExit { .. } => false,
        }
    }
}
