use compact_genome::interface::alphabet::AlphabetCharacter;
use strategies::{
    node_ord::NodeOrdStrategy, AlignmentStrategies, AlignmentStrategy, AlignmentStrategySelector,
};

use crate::{config::TemplateSwitchConfig, costs::cost::Cost};

use super::AlignmentGraphNode;

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
        /// The index that does not jump.
        primary_index: usize,
        /// The index that jumps.
        secondary_index: usize,
        gap_type: GapType,
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
        origin: TemplateSwitchPrimary,
        target: TemplateSwitchSecondary,
        first_offset: isize,
    },
    /// This node is the root node, hence it was not generated via alignment.
    Root,
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
                            context.offset1_costs.evaluate(&0),
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
                            context.offset1_costs.evaluate(&0),
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
                            context.offset1_costs.evaluate(&0),
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
                let (target_entrance_index, target_length) = match template_switch_secondary {
                    TemplateSwitchSecondary::Reference => {
                        (entrance_reference_index, reference.len())
                    }
                    TemplateSwitchSecondary::Query => (entrance_query_index, query.len()),
                };

                if template_switch_first_offset >= 0
                    && target_entrance_index as isize + template_switch_first_offset
                        < target_length as isize
                {
                    let old_cost = context
                        .offset1_costs
                        .evaluate(&template_switch_first_offset);
                    let new_cost = context
                        .offset1_costs
                        .evaluate(&(&template_switch_first_offset + 1));
                    assert!(new_cost >= old_cost);

                    let cost_increment = new_cost - old_cost;
                    output.extend(self.generate_template_switch_entrance_successor(
                        cost_increment,
                        template_switch_first_offset + 1,
                        context,
                    ))
                }

                if template_switch_first_offset <= 0
                    && target_entrance_index as isize + template_switch_first_offset > 0
                {
                    let old_cost = context
                        .offset1_costs
                        .evaluate(&template_switch_first_offset);
                    let new_cost = context
                        .offset1_costs
                        .evaluate(&(&template_switch_first_offset - 1));
                    assert!(new_cost >= old_cost);

                    let cost_increment = new_cost - old_cost;
                    output.extend(self.generate_template_switch_entrance_successor(
                        cost_increment,
                        template_switch_first_offset - 1,
                        context,
                    ))
                }

                output.extend(self.generate_secondary_root_node(context));
            }

            #[expect(unused)]
            Identifier::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                secondary_index,
                gap_type,
            } => {
                todo!()
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
                    template_switch_primary: template_switch_origin,
                    template_switch_secondary: template_switch_target,
                    template_switch_first_offset,
                    ..
                } => AlignmentType::TemplateSwitchEntrance {
                    origin: template_switch_origin,
                    target: template_switch_target,
                    first_offset: template_switch_first_offset,
                },

                Identifier::Secondary {
                    template_switch_primary,
                    template_switch_secondary,
                    primary_index,
                    secondary_index,
                    gap_type,
                    ..
                } => match gap_type {
                    GapType::Insertion => AlignmentType::Insertion,
                    GapType::Deletion => AlignmentType::Deletion,
                    GapType::None => {
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
                },
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
            template_switch_primary: template_switch_origin,
            template_switch_secondary: template_switch_target,
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
                    template_switch_primary: template_switch_origin,
                    template_switch_secondary: template_switch_target,
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
