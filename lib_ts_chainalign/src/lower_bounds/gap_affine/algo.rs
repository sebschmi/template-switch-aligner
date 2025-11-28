use std::fmt::Display;

use generic_a_star::{
    AStarContext, AStarNode,
    cost::{AStarCost, OrderedPairCost, U32Cost},
    reset::Reset,
};
use num_traits::Zero;

use crate::lower_bounds::gap_affine::GapAffineLowerBoundCostTable;

pub struct Context<'a, Cost> {
    costs: &'a GapAffineLowerBoundCostTable<Cost>,
    max_match_run: u32,
    max_n: usize,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    pub identifier: Identifier,
    pub cost: Cost,
    pub match_run: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Identifier {
    pub a: usize,
    pub b: usize,
    /// True if this node was reached via at least one non-match.
    pub has_non_match: bool,
    gap_type: GapType,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum GapType {
    None,
    InA,
    InB,
}

impl<'a, Cost> Context<'a, Cost> {
    pub fn new(
        costs: &'a GapAffineLowerBoundCostTable<Cost>,
        max_match_run: u32,
        max_n: usize,
    ) -> Self {
        Self {
            costs,
            max_match_run,
            max_n,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier {
                a: 0,
                b: 0,
                has_non_match: false,
                gap_type: GapType::None,
            },
            cost: Cost::zero(),
            match_run: 0,
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let Node {
            identifier:
                Identifier {
                    a,
                    b,
                    has_non_match,
                    gap_type,
                },
            cost,
            match_run,
        } = node;

        if *a < self.max_n && *b < self.max_n {
            if *match_run < self.max_match_run {
                // Match
                let new_cost = *cost;
                output.extend(std::iter::once(Node {
                    identifier: Identifier {
                        a: a + 1,
                        b: b + 1,
                        has_non_match: *has_non_match,
                        gap_type: GapType::None,
                    },
                    cost: new_cost,
                    match_run: match_run + 1,
                }));
            }

            // Substitution
            let new_cost = *cost + self.costs.substitution;
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    a: a + 1,
                    b: b + 1,
                    has_non_match: true,
                    gap_type: GapType::None,
                },
                cost: new_cost,
                match_run: 0,
            }));
        }

        if *a < self.max_n {
            // Gap in B
            let new_cost = *cost
                + match gap_type {
                    GapType::InB => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    a: a + 1,
                    b: *b,
                    has_non_match: true,
                    gap_type: GapType::InB,
                },
                cost: new_cost,
                match_run: 0,
            }));
        }

        if *b < self.max_n {
            // Gap in A
            let new_cost = *cost
                + match gap_type {
                    GapType::InA => self.costs.gap_extend,
                    _ => self.costs.gap_open,
                };
            output.extend(std::iter::once(Node {
                identifier: Identifier {
                    a: *a,
                    b: b + 1,
                    has_non_match: true,
                    gap_type: GapType::InA,
                },
                cost: new_cost,
                match_run: 0,
            }));
        }
    }

    fn is_target(&self, _node: &Self::Node) -> bool {
        // Run until whole matrix is filled
        false
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<Cost> Reset for Context<'_, Cost> {
    fn reset(&mut self) {
        // No internal state to reset
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = ();

    // Use match run as secondary cost
    type Cost = OrderedPairCost<Cost, U32Cost>;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Self::Cost {
        OrderedPairCost(self.cost, U32Cost::from_primitive(self.match_run))
    }

    fn a_star_lower_bound(&self) -> Self::Cost {
        OrderedPairCost(Cost::zero(), U32Cost::zero())
    }

    fn secondary_maximisable_score(&self) -> usize {
        0
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        // Backtracking not supported
        None
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        // Backtracking not supported
        None
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {}, {})", self.a, self.b, self.gap_type)
    }
}

impl Display for GapType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GapType::None => write!(f, "M/S"),
            GapType::InA => write!(f, "GA"),
            GapType::InB => write!(f, "GB"),
        }
    }
}

impl<Cost: Ord> PartialOrd for Node<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Cost: Ord> Ord for Node<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost
            .cmp(&other.cost)
            .then_with(|| self.match_run.cmp(&other.match_run))
    }
}
