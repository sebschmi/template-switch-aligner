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
        template_switch_origin: TemplateSwitchOrigin,
        template_switch_target: TemplateSwitchTarget,
        template_switch_first_offset: isize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GapType {
    Insertion,
    Deletion,
    None,
}

/// Origin is the sequence that gets replaced.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TemplateSwitchOrigin {
    Reference,
    Query,
}

/// Target is the sequence where the replacement comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TemplateSwitchTarget {
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
        origin: TemplateSwitchOrigin,
        target: TemplateSwitchTarget,
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
                        output.extend([self.generate_primary_diagonal_successor(
                            0,
                            context
                                .primary_edit_costs
                                .match_or_substitution_cost(r.clone(), q.clone()),
                            context,
                        )]);
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend([self.generate_primary_diagonal_successor(
                            flank_index + 1,
                            context
                                .left_flank_edit_costs
                                .match_or_substitution_cost(r, q),
                            context,
                        )]);
                    } else if flank_index < 0 {
                        output.extend([self.generate_primary_diagonal_successor(
                            flank_index + 1,
                            context
                                .right_flank_edit_costs
                                .match_or_substitution_cost(r, q),
                            context,
                        )]);
                    } else if flank_index == context.left_flank_length {
                        todo!("generate template switch start nodes");
                    }
                }

                if reference_index < reference.len() {
                    // Deleted character
                    let r = reference[reference_index].clone();

                    if flank_index == 0 {
                        output.extend([self.generate_primary_deletion_successor(
                            0,
                            context
                                .primary_edit_costs
                                .gap_costs(r.clone(), gap_type != GapType::Deletion),
                            context,
                        )]);
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend([self.generate_primary_deletion_successor(
                            flank_index + 1,
                            context
                                .left_flank_edit_costs
                                .gap_costs(r, gap_type != GapType::Deletion),
                            context,
                        )]);
                    } else if flank_index < 0 {
                        output.extend([self.generate_primary_deletion_successor(
                            flank_index + 1,
                            context
                                .right_flank_edit_costs
                                .gap_costs(r, gap_type != GapType::Deletion),
                            context,
                        )]);
                    } else if flank_index == context.left_flank_length {
                        todo!("generate template switch start nodes");
                    }
                }

                if query_index < query.len() {
                    // Inserted character
                    let q = query[query_index].clone();

                    if flank_index == 0 {
                        output.extend([self.generate_primary_insertion_successor(
                            0,
                            context
                                .primary_edit_costs
                                .gap_costs(q.clone(), gap_type != GapType::Insertion),
                            context,
                        )]);
                    }

                    if flank_index >= 0 && flank_index < context.left_flank_length {
                        output.extend([self.generate_primary_insertion_successor(
                            flank_index + 1,
                            context
                                .left_flank_edit_costs
                                .gap_costs(q, gap_type != GapType::Insertion),
                            context,
                        )]);
                    } else if flank_index < 0 {
                        output.extend([self.generate_primary_insertion_successor(
                            flank_index + 1,
                            context
                                .right_flank_edit_costs
                                .gap_costs(q, gap_type != GapType::Insertion),
                            context,
                        )]);
                    } else if flank_index == context.left_flank_length {
                        todo!("generate template switch start nodes");
                    }
                }
            }

            #[expect(unused)]
            Identifier::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_origin,
                template_switch_target,
                template_switch_first_offset,
            } => {
                let (origin_entrance_index, origin_length) = match template_switch_origin {
                    TemplateSwitchOrigin::Reference => (entrance_reference_index, reference.len()),
                    TemplateSwitchOrigin::Query => (entrance_query_index, query.len()),
                };
                let (target_entrance_index, target_length) = match template_switch_target {
                    TemplateSwitchTarget::Reference => (entrance_reference_index, reference.len()),
                    TemplateSwitchTarget::Query => (entrance_query_index, query.len()),
                };

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
                    template_switch_origin,
                    template_switch_target,
                    template_switch_first_offset,
                    ..
                } => AlignmentType::TemplateSwitchEntrance {
                    origin: template_switch_origin,
                    target: template_switch_target,
                    first_offset: template_switch_first_offset,
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
    ) -> Self {
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
        Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        }
    }

    fn generate_primary_deletion_successor(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Self {
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
        Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        }
    }

    fn generate_primary_insertion_successor(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &<Self as AlignmentGraphNode<Strategies::Alphabet>>::Context,
    ) -> Self {
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
        Self {
            node_data,
            strategies: self.strategies.generate_successor(context),
        }
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
