use std::fmt::Display;

use generic_a_star::{AStarContext, AStarNode, cost::AStarCost, reset::Reset};

use crate::{
    alignment::{
        AlignmentType, GapType, coordinates::AlignmentCoordinates, sequences::AlignmentSequences,
    },
    costs::GapAffineCosts,
};

pub struct Context<'costs, 'sequences, Cost> {
    costs: &'costs GapAffineCosts<Cost>,
    sequences: &'sequences AlignmentSequences,
    start: AlignmentCoordinates,
    end: AlignmentCoordinates,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub predecessor: Option<Identifier>,
    pub predecessor_alignment_type: Option<AlignmentType>,
    pub cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Identifier {
    pub coordinates: AlignmentCoordinates,
    gap_type: GapType,
}

impl<'costs, 'sequences, Cost> Context<'costs, 'sequences, Cost> {
    pub fn new(
        costs: &'costs GapAffineCosts<Cost>,
        sequences: &'sequences AlignmentSequences,
        start: AlignmentCoordinates,
        end: AlignmentCoordinates,
    ) -> Self {
        Self {
            costs,
            sequences,
            start,
            end,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                coordinates: self.start,
                gap_type: GapType::None,
            },
            predecessor: None,
            predecessor_alignment_type: None,
            cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let Node {
            identifier, cost, ..
        } = node;
        let predecessor = Some(*identifier);
        let Identifier {
            coordinates,
            gap_type,
        } = *identifier;

        if coordinates.can_increment_both(self.start, self.end) {
            let (c1, c2) = self.sequences.characters(coordinates);
            let is_match = c1 == c2;

            if is_match {
                // Match
                let new_cost = *cost;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        gap_type: GapType::None,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::Match),
                    cost: new_cost,
                }));
            } else {
                // Substitution
                let new_cost = *cost + self.costs.substitution;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        coordinates: coordinates.increment_both(),
                        gap_type: GapType::None,
                    },
                    predecessor,
                    predecessor_alignment_type: Some(AlignmentType::Substitution),
                    cost: new_cost,
                }));
            }
        }

        if coordinates.can_increment_1(self.start, self.end) {
            // Gap in 2
            let new_cost = *cost
                + match gap_type {
                    GapType::In2 => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_1(),
                    gap_type: GapType::In2,
                },
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::Gap2),
                cost: new_cost,
            }));
        }

        if coordinates.can_increment_2(self.start, self.end) {
            // Gap in 1
            let new_cost = *cost
                + match gap_type {
                    GapType::In1 => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    coordinates: coordinates.increment_2(),
                    gap_type: GapType::In1,
                },
                predecessor,
                predecessor_alignment_type: Some(AlignmentType::Gap1),
                cost: new_cost,
            }));
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier.coordinates == self.end
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<Cost> Reset for Context<'_, '_, Cost> {
    fn reset(&mut self) {
        unimplemented!()
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = AlignmentType;

    // Use match run as secondary cost
    type Cost = Cost;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        Cost::zero()
    }

    fn secondary_maximisable_score(&self) -> usize {
        0
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        self.predecessor_alignment_type
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {})",
            self.coordinates.seq1(),
            self.coordinates.seq2(),
            self.gap_type
        )
    }
}

impl<Cost: Ord> PartialOrd for Node<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Cost: Ord> Ord for Node<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.cmp(&other.cost)
    }
}
