use std::marker::PhantomData;

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{reset::Reset, AStarContext, AStarNode};

use crate::costs::cost::Cost;

use super::{alignment_result::IAlignmentType, AlignmentContext};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Node {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    predecessor_edge_type: AlignmentType,
    cost: Cost,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Identifier {
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

#[derive(Debug, Clone)]
pub struct Context<
    'reference,
    'query,
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
> {
    reference: &'reference SubsequenceType,
    query: &'query SubsequenceType,

    scoring_table: ScoringTable,
    phantom_data: PhantomData<AlphabetType>,
}

impl AStarNode for Node {
    type Identifier = Identifier;

    type EdgeType = AlignmentType;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Cost {
        0.into()
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        Some(self.predecessor_edge_type)
    }
}

impl<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    > AStarContext for Context<'_, '_, AlphabetType, SubsequenceType>
{
    type Node = Node;

    fn create_root(&self) -> Self::Node {
        Self::Node {
            identifier: Identifier::new(0, 0, GapType::None),
            predecessor: None,
            predecessor_edge_type: AlignmentType::Root,
            cost: Cost::ZERO,
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        if node.identifier.reference_index < self.reference.len()
            && node.identifier.query_index < self.query.len()
        {
            let is_match = self.reference[node.identifier.reference_index]
                == self.query[node.identifier.query_index];
            output.extend([Self::Node {
                identifier: node.identifier.increment_both(),
                predecessor: Some(node.identifier),
                predecessor_edge_type: if is_match {
                    AlignmentType::Match
                } else {
                    AlignmentType::Substitution
                },
                cost: node.cost
                    + if is_match {
                        self.scoring_table.match_cost
                    } else {
                        self.scoring_table.substitution_cost
                    },
            }]);
        }

        if node.identifier.reference_index < self.reference.len() {
            output.extend([Self::Node {
                identifier: node.identifier.increment_reference(),
                predecessor: Some(node.identifier),
                predecessor_edge_type: AlignmentType::Deletion,
                cost: node.cost
                    + if let Some(predecessor) = node.predecessor {
                        if predecessor.increment_reference() == node.identifier {
                            self.scoring_table.gap_extend_cost
                        } else {
                            self.scoring_table.gap_open_cost
                        }
                    } else {
                        self.scoring_table.gap_open_cost
                    },
            }]);
        }

        if node.identifier.query_index < self.query.len() {
            output.extend([Self::Node {
                identifier: node.identifier.increment_query(),
                predecessor: Some(node.identifier),
                predecessor_edge_type: AlignmentType::Insertion,
                cost: node.cost
                    + if let Some(predecessor) = node.predecessor {
                        if predecessor.increment_query() == node.identifier {
                            self.scoring_table.gap_extend_cost
                        } else {
                            self.scoring_table.gap_open_cost
                        }
                    } else {
                        self.scoring_table.gap_open_cost
                    },
            }]);
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.reference_index == self.reference.len()
            && node.identifier.query_index == self.query.len()
    }

    fn max_cost(&self) -> Option<Cost> {
        None
    }
}

impl<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    > Reset for Context<'_, '_, AlphabetType, SubsequenceType>
{
    fn reset(&mut self) {}
}

impl<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    > AlignmentContext for Context<'_, '_, AlphabetType, SubsequenceType>
{
    type AlphabetType = AlphabetType;

    type AlignmentType = AlignmentType;

    type SubsequenceType = SubsequenceType;

    fn reference(&self) -> &Self::SubsequenceType {
        self.reference
    }

    fn query(&self) -> &Self::SubsequenceType {
        self.query
    }
}

impl<
        'reference,
        'query,
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    > Context<'reference, 'query, AlphabetType, SubsequenceType>
{
    pub fn new(
        reference: &'reference SubsequenceType,
        query: &'query SubsequenceType,
        scoring_table: ScoringTable,
    ) -> Self {
        Self {
            reference,
            query,
            scoring_table,
            phantom_data: PhantomData,
        }
    }
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
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

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            identifier,
            predecessor,
            predecessor_edge_type,
            cost,
        } = self;
        write!(f, "{identifier}; ")?;
        if let Some(predecessor) = predecessor {
            write!(f, "predecessor: {predecessor}; ")?;
        }
        write!(f, "alignment_type: {predecessor_edge_type}")?;
        write!(f, "cost: {cost}")
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

    const fn anti_diagonal(&self) -> usize {
        self.reference_index + self.query_index
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
        write!(
            f,
            "({}, {}, {})",
            self.reference_index, self.query_index, self.gap_type
        )
    }
}

impl IAlignmentType for AlignmentType {
    fn is_repeatable(&self) -> bool {
        true
    }

    fn is_repeated(&self, previous: &Self) -> bool {
        self == previous
    }

    fn is_internal(&self) -> bool {
        self == &Self::Root
    }
}
