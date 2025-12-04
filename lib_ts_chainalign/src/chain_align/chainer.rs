use std::{fmt::Display, iter};

use generic_a_star::{AStarContext, AStarNode, cost::AStarCost, reset::Reset};

use crate::{anchors::Anchors, chaining_cost_function::ChainingCostFunction};

const DEBUG_CHAINER: bool = false;

pub struct Context<'anchors, 'chaining_cost_function, Cost> {
    anchors: &'anchors Anchors,
    chaining_cost_function: &'chaining_cost_function mut ChainingCostFunction<Cost>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Node<Cost> {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum Identifier {
    Start,
    Primary { index: usize },
    End,
}

impl<'anchors, 'chaining_cost_function, Cost> Context<'anchors, 'chaining_cost_function, Cost> {
    pub fn new(
        anchors: &'anchors Anchors,
        chaining_cost_function: &'chaining_cost_function mut ChainingCostFunction<Cost>,
    ) -> Self {
        Self {
            anchors,
            chaining_cost_function,
        }
    }
}

impl<Cost: AStarCost> AStarContext for Context<'_, '_, Cost> {
    type Node = Node<Cost>;

    fn create_root(&self) -> Self::Node {
        Node {
            identifier: Identifier::Start,
            predecessor: None,
            cost: Cost::zero(),
        }
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        let predecessor = Some(node.identifier);
        let predecessor_cost = node.cost;

        if DEBUG_CHAINER {
            println!("Generating successors of {node}");
        }

        match node.identifier {
            Identifier::Start => {
                output.extend(
                    (0..self.anchors.primary.len())
                        .flat_map(|successor_index| {
                            if DEBUG_CHAINER {
                                println!(
                                    "Checking anchor P-{successor_index}: {}",
                                    self.anchors.primary[successor_index]
                                );
                            }

                            let cost = predecessor_cost.checked_add(
                                &self
                                    .chaining_cost_function
                                    .primary_from_start(successor_index),
                            )?;
                            if DEBUG_CHAINER {
                                println!("Cost: {cost}");
                            }

                            (cost != Cost::max_value()).then_some(Node {
                                identifier: Identifier::Primary {
                                    index: successor_index,
                                },
                                predecessor,
                                cost,
                            })
                        })
                        .chain(iter::once({
                            let cost = predecessor_cost
                                .checked_add(&self.chaining_cost_function.start_to_end())
                                .unwrap();
                            debug_assert_ne!(cost, Cost::max_value());
                            Node {
                                identifier: Identifier::End,
                                predecessor,
                                cost,
                            }
                        })),
                );
            }
            Identifier::Primary { index } => output.extend(
                (0..self.anchors.primary.len())
                    .flat_map(|successor_index| {
                        if DEBUG_CHAINER {
                            println!(
                                "Checking anchor P-{successor_index}: {}",
                                self.anchors.primary[successor_index]
                            );
                        }

                        let cost = predecessor_cost.checked_add(
                            &self.chaining_cost_function.primary(index, successor_index),
                        )?;
                        if DEBUG_CHAINER {
                            println!("Cost: {cost}");
                        }

                        (cost != Cost::max_value()).then_some(Node {
                            identifier: Identifier::Primary {
                                index: successor_index,
                            },
                            predecessor,
                            cost,
                        })
                    })
                    .chain(iter::once({
                        let cost = predecessor_cost
                            .checked_add(&self.chaining_cost_function.primary_to_end(index))
                            .unwrap();
                        debug_assert_ne!(cost, Cost::max_value());
                        Node {
                            identifier: Identifier::End,
                            predecessor,
                            cost,
                        }
                    })),
            ),
            Identifier::End => { /* Has no successors */ }
        }
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        node.identifier == Identifier::End
    }

    fn cost_limit(&self) -> Option<<Self::Node as generic_a_star::AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<Cost: AStarCost> Reset for Context<'_, '_, Cost> {
    fn reset(&mut self) {
        // Nothing to do.
    }
}

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = Identifier;

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
        self.predecessor
    }
}

impl<Cost: Ord> Ord for Node<Cost> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.cmp(&other.cost)
    }
}

impl<Cost: Ord> PartialOrd for Node<Cost> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Cost: Display> Display for Node<Cost> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.identifier, self.cost)
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Identifier::Start => write!(f, "start"),
            Identifier::Primary { index } => write!(f, "P-{index}"),
            Identifier::End => write!(f, "end"),
        }
    }
}
