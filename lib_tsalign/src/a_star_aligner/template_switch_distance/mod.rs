use crate::cost::Cost;

use super::AlignmentGraphNode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Node {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Identifier {
    Primary {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum GapType {
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

impl AlignmentGraphNode for Node {
    type Identifier = Identifier;

    type Context = ScoringTable;

    type AlignmentType = AlignmentType;

    fn create_root() -> Self {
        Self {
            identifier: Identifier::new_primary(0, 0, GapType::None),
            predecessor: None,
            cost: Cost::ZERO,
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
        match self.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                gap_type,
            } => {
                if reference_index < reference.len() && query_index < query.len() {
                    output.extend([Self {
                        identifier: self.identifier.increment_both(),
                        predecessor: Some(self.identifier),
                        cost: self.cost
                            + if reference[reference_index] == query[query_index] {
                                context.match_cost
                            } else {
                                context.substitution_cost
                            },
                    }]);
                }

                if reference_index < reference.len() {
                    let identifier = self.identifier.increment_reference();
                    output.extend([Self {
                        identifier,
                        predecessor: Some(self.identifier),
                        cost: self.cost
                            + if gap_type == identifier.gap_type() {
                                context.gap_extend_cost
                            } else {
                                context.gap_open_cost
                            },
                    }]);
                }

                if query_index < query.len() {
                    let identifier = self.identifier.increment_query();
                    output.extend([Self {
                        identifier,
                        predecessor: Some(self.identifier),
                        cost: self.cost
                            + if gap_type == identifier.gap_type() {
                                context.gap_extend_cost
                            } else {
                                context.gap_open_cost
                            },
                    }]);
                }
            }
        }
    }

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Cost {
        self.cost
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
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
        if self.predecessor.is_some() {
            match self.identifier {
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
        match self.identifier {
            Identifier::Primary {
                reference_index,
                query_index,
                ..
            } => reference_index == reference.len() && query_index == query.len(),
        }
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.cost.cmp(&other.cost) {
            // This secondary ordering may make things actually slower.
            // While it does reduce the number of visited nodes a little bit,
            // it also makes heap operations more expensive.
            // Preliminary testing showed that this would be a slowdown.
            std::cmp::Ordering::Equal => other
                .identifier
                .anti_diagonal()
                .cmp(&self.identifier.anti_diagonal()),
            ordering => ordering,
        }
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

impl PartialOrd for Node {
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
