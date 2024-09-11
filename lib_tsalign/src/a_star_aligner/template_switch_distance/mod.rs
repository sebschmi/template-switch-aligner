use strategies::{
    node_ord::NodeOrdStrategy, AlignmentStrategies, AlignmentStrategy, AlignmentStrategySelector,
};

use crate::{cost::Cost, cost_table::TemplateSwitchCostTable};

use super::AlignmentGraphNode;

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
        steps_since_last_template_switch: usize,
        gap_type: GapType,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GapType {
    Insertion,
    Deletion,
    None,
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
    /// This node is the root node, hence it was not generated via alignment.
    Root,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Context<Alphabet> {
    pub costs: TemplateSwitchCostTable<Alphabet>,
    pub left_flank_length: usize,
    pub right_flank_length: usize,
}

impl<Strategies: AlignmentStrategySelector> AlignmentGraphNode<Strategies::Alphabet>
    for Node<Strategies>
{
    type Identifier = Identifier;

    type Context = Context<Strategies::Alphabet>;

    type AlignmentType = AlignmentType;

    fn create_root(context: &Self::Context) -> Self {
        Self {
            node_data: NodeData {
                identifier: Identifier::new_primary(0, 0, usize::MAX, GapType::None),
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
                gap_type,
                ..
            } => {
                if reference_index < reference.len() && query_index < query.len() {
                    output.extend([Self {
                        node_data: NodeData {
                            identifier: self.node_data.identifier.generate_diagonal_successor(),
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + context.costs.primary_edit_costs.match_or_substitution_cost(
                                    reference[reference_index].clone(),
                                    query[query_index].clone(),
                                ),
                        },
                        strategies: self.strategies.generate_successor(context),
                    }]);
                }

                if reference_index < reference.len() {
                    let identifier = self.node_data.identifier.generate_deletion_successor();
                    output.extend([Self {
                        node_data: NodeData {
                            identifier,
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + if gap_type == identifier.gap_type() {
                                    context
                                        .costs
                                        .primary_edit_costs
                                        .gap_extend_cost(reference[reference_index].clone())
                                } else {
                                    context
                                        .costs
                                        .primary_edit_costs
                                        .gap_open_cost(reference[reference_index].clone())
                                },
                        },
                        strategies: self.strategies.generate_successor(context),
                    }]);
                }

                if query_index < query.len() {
                    let identifier = self.node_data.identifier.generate_insertion_successor();
                    output.extend([Self {
                        node_data: NodeData {
                            identifier,
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + if gap_type == identifier.gap_type() {
                                    context
                                        .costs
                                        .primary_edit_costs
                                        .gap_extend_cost(query[query_index].clone())
                                } else {
                                    context
                                        .costs
                                        .primary_edit_costs
                                        .gap_open_cost(query[query_index].clone())
                                },
                        },
                        strategies: self.strategies.generate_successor(context),
                    }]);
                }
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

impl Identifier {
    const fn new_primary(
        reference_index: usize,
        query_index: usize,
        steps_since_last_template_switch: usize,
        gap_type: GapType,
    ) -> Self {
        Self::Primary {
            reference_index,
            query_index,
            steps_since_last_template_switch,
            gap_type,
        }
    }

    const fn generate_deletion_successor(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                steps_since_last_template_switch,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index,
                steps_since_last_template_switch: steps_since_last_template_switch
                    .saturating_add(1),
                gap_type: GapType::Deletion,
            },
        }
    }

    const fn generate_insertion_successor(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                steps_since_last_template_switch,
                ..
            } => Self::Primary {
                reference_index,
                query_index: query_index + 1,
                steps_since_last_template_switch: steps_since_last_template_switch
                    .saturating_add(1),
                gap_type: GapType::Insertion,
            },
        }
    }

    const fn generate_diagonal_successor(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                steps_since_last_template_switch,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index: query_index + 1,
                steps_since_last_template_switch: steps_since_last_template_switch
                    .saturating_add(1),
                gap_type: GapType::None,
            },
        }
    }

    const fn gap_type(self) -> GapType {
        match self {
            Self::Primary { gap_type, .. } => gap_type,
        }
    }

    /// Returns the anti-diagonal for variants where it exists, or zero otherwise.
    const fn anti_diagonal(self) -> usize {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => reference_index + query_index,
        }
    }
}

impl<Strategies: AlignmentStrategySelector> PartialOrd for Node<Strategies> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::fmt::Display for AlignmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentType::Insertion => write!(f, "I"),
            AlignmentType::Deletion => write!(f, "D"),
            AlignmentType::Substitution => write!(f, "S"),
            AlignmentType::Match => write!(f, "M"),
            AlignmentType::Root => Ok(()),
        }
    }
}

impl std::fmt::Display for GapType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GapType::Insertion => write!(f, "I"),
            GapType::Deletion => write!(f, "D"),
            GapType::None => write!(f, "M/S"),
        }
    }
}

impl std::fmt::Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                steps_since_last_template_switch,
                gap_type,
            } => write!(
                f,
                "Primary({}, {}, {}, {})",
                reference_index, query_index, steps_since_last_template_switch, gap_type
            ),
        }
    }
}
