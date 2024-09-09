use strategies::{
    node_ord::NodeOrdStrategy, AlignmentStrategies, AlignmentStrategy, AlignmentStrategySelector,
};

use crate::cost::Cost;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScoringTable {
    pub match_cost: Cost,
    pub substitution_cost: Cost,
    pub gap_open_cost: Cost,
    pub gap_extend_cost: Cost,
}

impl<Strategies: AlignmentStrategySelector> AlignmentGraphNode for Node<Strategies> {
    type Identifier = Identifier;

    type Context = ScoringTable;

    type AlignmentType = AlignmentType;

    fn create_root() -> Self {
        Self {
            node_data: NodeData {
                identifier: Identifier::new_primary(0, 0, GapType::None),
                predecessor: None,
                cost: Cost::ZERO,
            },
            strategies: AlignmentStrategy::create_root(),
        }
    }

    fn generate_successors<
        AlphabetType: compact_genome::interface::alphabet::Alphabet,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
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
            } => {
                if reference_index < reference.len() && query_index < query.len() {
                    output.extend([Self {
                        node_data: NodeData {
                            identifier: self.node_data.identifier.increment_both(),
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + if reference[reference_index] == query[query_index] {
                                    context.match_cost
                                } else {
                                    context.substitution_cost
                                },
                        },
                        strategies: self.strategies.generate_successor(),
                    }]);
                }

                if reference_index < reference.len() {
                    let identifier = self.node_data.identifier.increment_reference();
                    output.extend([Self {
                        node_data: NodeData {
                            identifier,
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + if gap_type == identifier.gap_type() {
                                    context.gap_extend_cost
                                } else {
                                    context.gap_open_cost
                                },
                        },
                        strategies: self.strategies.generate_successor(),
                    }]);
                }

                if query_index < query.len() {
                    let identifier = self.node_data.identifier.increment_query();
                    output.extend([Self {
                        node_data: NodeData {
                            identifier,
                            predecessor: Some(self.node_data.identifier),
                            cost: self.node_data.cost
                                + if gap_type == identifier.gap_type() {
                                    context.gap_extend_cost
                                } else {
                                    context.gap_open_cost
                                },
                        },
                        strategies: self.strategies.generate_successor(),
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
        AlphabetType: compact_genome::interface::alphabet::Alphabet,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
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
        AlphabetType: compact_genome::interface::alphabet::Alphabet,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
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
    const fn new_primary(reference_index: usize, query_index: usize, gap_type: GapType) -> Self {
        Self::Primary {
            reference_index,
            query_index,
            gap_type,
        }
    }

    const fn increment_reference(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index,
                gap_type: GapType::Deletion,
            },
        }
    }

    const fn increment_query(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index,
                query_index: query_index + 1,
                gap_type: GapType::Insertion,
            },
        }
    }

    const fn increment_both(self) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => Self::Primary {
                reference_index: reference_index + 1,
                query_index: query_index + 1,
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
                gap_type,
            } => write!(
                f,
                "Primary({}, {}, {})",
                reference_index, query_index, gap_type
            ),
        }
    }
}
