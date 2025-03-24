use generic_a_star::{AStarContext, AStarNode, cost::AStarCost, reset::Reset};

use crate::seed::ChainingAnchors;

use super::node::{Identifier, Node};

pub trait ChainingCostsProvider {
    type Cost: AStarCost;

    /// Returns the cost of chaining `to` after `from`.
    ///
    /// If the chaining is impossible, then `Cost::MAX` is returned.
    fn chaining_costs(&self, from: &Identifier, to: &Identifier) -> Self::Cost;
}

pub struct Context<ChainingCosts: ChainingCostsProvider> {
    chaining_costs: ChainingCosts,
    chaining_anchors: ChainingAnchors,
}

impl<ChainingCosts: ChainingCostsProvider> AStarContext for Context<ChainingCosts> {
    type Node = Node<<ChainingCosts as ChainingCostsProvider>::Cost>;

    fn create_root(&self) -> Self::Node {
        Self::Node::new_root()
    }

    fn generate_successors(&mut self, node: &Self::Node, output: &mut impl Extend<Self::Node>) {
        // We skip all anchors whose reference block is not right of the current node.
        let first_chainable_block_index = match node.identifier() {
            Identifier::Root => 0,
            Identifier::Anchor { anchor } => {
                self.chaining_anchors
                    .anchors()
                    .partition_point(|chaining_anchor| {
                        // We do not allow to chain overlapping blocks.
                        chaining_anchor.reference_block().start < anchor.reference_block().end
                    })
            }
            Identifier::Target => self.chaining_anchors.anchors().len(),
        };

        output.extend(
            self.chaining_anchors
                .anchors()
                .iter()
                .skip(first_chainable_block_index)
                .filter_map(|chaining_anchor| {
                    // We skip all anchors whose query block is not right of the current node.
                    if let Identifier::Anchor { anchor } = node.identifier() {
                        if chaining_anchor.query_block().start < anchor.query_block().end {
                            return None;
                        }
                    }

                    let successor_identifier = Identifier::Anchor {
                        anchor: chaining_anchor.clone(),
                    };
                    let cost_increment = self
                        .chaining_costs
                        .chaining_costs(node.identifier(), &successor_identifier);
                    node.generate_successor(successor_identifier, cost_increment)
                })
                .chain(if !matches!(node.identifier(), Identifier::Target) {
                    node.generate_successor(
                        Identifier::Target,
                        self.chaining_costs
                            .chaining_costs(node.identifier(), &Identifier::Target),
                    )
                } else {
                    None
                }),
        );
    }

    fn is_target(&self, node: &Self::Node) -> bool {
        matches!(node.identifier(), Identifier::Target)
    }

    fn cost_limit(&self) -> Option<<Self::Node as AStarNode>::Cost> {
        None
    }

    fn memory_limit(&self) -> Option<usize> {
        None
    }
}

impl<ChainingCosts: ChainingCostsProvider> Context<ChainingCosts> {
    pub fn new(chaining_costs: ChainingCosts, chaining_anchors: ChainingAnchors) -> Self {
        Self {
            chaining_costs,
            chaining_anchors,
        }
    }
}

impl<ChainingCosts: ChainingCostsProvider> Reset for Context<ChainingCosts> {
    fn reset(&mut self) {
        todo!()
    }
}
