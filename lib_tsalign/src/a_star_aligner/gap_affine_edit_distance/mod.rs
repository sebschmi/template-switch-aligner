use std::fmt::Display;

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

use crate::cost::Cost;

use super::AlignmentGraphNode;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Node {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct Identifier {
    reference_index: usize,
    query_index: usize,
    gap_type: GapType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum GapType {
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
            identifier: Identifier::new(0, 0, GapType::None),
            predecessor: None,
            cost: 0.into(),
        }
    }

    fn generate_successors<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        context: &Self::Context,
        output: &mut impl Extend<Self>,
    ) {
        if self.identifier.reference_index < reference.len()
            && self.identifier.query_index < query.len()
        {
            output.extend([Self {
                identifier: self.identifier.increment_both(),
                predecessor: Some(self.identifier),
                cost: self.cost
                    + if reference[self.identifier.reference_index]
                        == query[self.identifier.query_index]
                    {
                        context.match_cost
                    } else {
                        context.substitution_cost
                    },
            }]);
        }

        if self.identifier.reference_index < reference.len() {
            output.extend([Self {
                identifier: self.identifier.increment_reference(),
                predecessor: Some(self.identifier),
                cost: self.cost
                    + if let Some(predecessor) = self.predecessor {
                        if predecessor.increment_reference() == self.identifier {
                            context.gap_extend_cost
                        } else {
                            context.gap_open_cost
                        }
                    } else {
                        context.gap_open_cost
                    },
            }]);
        }

        if self.identifier.query_index < query.len() {
            output.extend([Self {
                identifier: self.identifier.increment_query(),
                predecessor: Some(self.identifier),
                cost: self.cost
                    + if let Some(predecessor) = self.predecessor {
                        if predecessor.increment_query() == self.identifier {
                            context.gap_extend_cost
                        } else {
                            context.gap_open_cost
                        }
                    } else {
                        context.gap_open_cost
                    },
            }]);
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
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        _context: &Self::Context,
    ) -> Self::AlignmentType {
        if let Some(predecessor) = self.predecessor {
            if predecessor.increment_both() == self.identifier {
                if reference[self.identifier.reference_index - 1]
                    == query[self.identifier.query_index - 1]
                {
                    AlignmentType::Match
                } else {
                    AlignmentType::Substitution
                }
            } else if predecessor.increment_reference() == self.identifier {
                AlignmentType::Deletion
            } else if predecessor.increment_query() == self.identifier {
                AlignmentType::Insertion
            } else {
                unreachable!("this node was generated from the predecessor in one of the three ways handled above")
            }
        } else {
            AlignmentType::Root
        }
    }

    fn is_target<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        reference: &SubsequenceType,
        query: &SubsequenceType,
        _context: &Self::Context,
    ) -> bool {
        self.identifier.reference_index == reference.len()
            && self.identifier.query_index == query.len()
    }
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        /*match self.cost.cmp(&other.cost) {
            std::cmp::Ordering::Equal => other
                .identifier
                .anti_diagonal()
                .cmp(&self.identifier.anti_diagonal()),
            ordering => ordering,
        }*/
        /*match self.cost.cmp(&other.cost) {
            std::cmp::Ordering::Equal => self
                .predecessor_gap_secondary_cost()
                .cmp(&other.predecessor_gap_secondary_cost()),
            ordering => ordering,
        }*/
        self.cost.cmp(&other.cost)
    }
}

impl Identifier {
    const fn new(reference_index: usize, query_index: usize, gap_type: GapType) -> Self {
        Self {
            reference_index,
            query_index,
            gap_type,
        }
    }

    const fn increment_reference(&self) -> Self {
        Self {
            reference_index: self.reference_index + 1,
            query_index: self.query_index,
            gap_type: GapType::Deletion,
        }
    }

    const fn increment_query(&self) -> Self {
        Self {
            reference_index: self.reference_index,
            query_index: self.query_index + 1,
            gap_type: GapType::Insertion,
        }
    }

    const fn increment_both(&self) -> Self {
        Self {
            reference_index: self.reference_index + 1,
            query_index: self.query_index + 1,
            gap_type: GapType::None,
        }
    }
}

impl Display for AlignmentType {
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

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.reference_index, self.query_index)
    }
}
