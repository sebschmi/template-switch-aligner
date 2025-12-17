use std::cmp::Ordering;

use compare::Compare;
use rustc_hash::{FxHashMapSeed, FxSeededState};

use crate::{AStarNode, comparator::AStarNodeComparator, reset::Reset};

/// The closed list for the A* algorithm.
pub trait AStarClosedList<Node: AStarNode>: Reset {
    /// Create a new empty closed list.
    fn new() -> Self;

    /// Returns the number of nodes in the closed list.
    fn len(&self) -> usize;

    /// Returns true if the closed list is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Insert a node with the given identifier into the closed list.
    ///
    /// If there is no node mapped to the given identifier, then `None` is returned.
    /// Otherwise, if the currently mapped node is better than the inserted node, then the inserted node is returned.
    /// And if the currently mapped node is worse than the inserted node, then the currently mapped node is replaced and returned.
    fn insert(&mut self, identifier: <Node as AStarNode>::Identifier, node: Node) -> Option<Node>;

    /// Return a reference to the node specified by `identifier`.
    ///
    /// Returns `None` if no node with the given identifier exists.
    fn get(&self, identifier: &<Node as AStarNode>::Identifier) -> Option<&Node>;

    /// Returns true if the given identifier is mapped to a node.
    fn contains_identifier(&self, identifier: &<Node as AStarNode>::Identifier) -> bool {
        self.get(identifier).is_some()
    }

    /// Returns an iterator over the nodes in the closed list.
    fn iter<'this: 'node, 'node>(
        &'this self,
    ) -> impl use<'this, 'node, Self, Node> + Iterator<Item = &'node Node>
    where
        Node: 'node;

    fn can_skip_node(&self, node: &Node, is_label_setting: bool) -> bool {
        if let Some(previous_visit) = self.get(node.identifier()) {
            if is_label_setting {
                // In label-setting mode, if we have already visited the node, we now must be visiting it with a higher or equal cost.
                debug_assert!(
                    previous_visit.cost() + previous_visit.a_star_lower_bound()
                        <= node.cost() + node.a_star_lower_bound(),
                    "Revisiting node {} (previous: {}) at lower costs:\n{}",
                    node,
                    previous_visit,
                    {
                        use std::fmt::Write;
                        let mut previous_visit = previous_visit;
                        let mut node = node;
                        let mut out = String::new();

                        writeln!(out, "previous_visit:").unwrap();
                        while let Some(predecessor) = previous_visit.predecessor() {
                            writeln!(out, "{previous_visit}").unwrap();
                            previous_visit = self.get(predecessor).unwrap();
                        }

                        writeln!(out, "\nnode:").unwrap();
                        while let Some(predecessor) = node.predecessor() {
                            writeln!(out, "{node}").unwrap();
                            node = self.get(predecessor).unwrap();
                        }

                        out
                    }
                );

                true
            } else if AStarNodeComparator.compare(node, previous_visit) != Ordering::Greater {
                // If we are label-correcting, we may still find a better node later on.
                // Skip if equal or worse.
                true
            } else {
                false
            }
        } else {
            false
        }
    }
}

impl<Node: AStarNode> AStarClosedList<Node>
    for FxHashMapSeed<<Node as AStarNode>::Identifier, Node>
{
    fn new() -> Self {
        FxHashMapSeed::with_hasher(FxSeededState::with_seed(0))
    }

    fn len(&self) -> usize {
        Self::len(self)
    }

    fn insert(&mut self, identifier: <Node as AStarNode>::Identifier, node: Node) -> Option<Node> {
        Self::insert(self, identifier, node)
    }

    fn get(&self, identifier: &<Node as AStarNode>::Identifier) -> Option<&Node> {
        Self::get(self, identifier)
    }

    fn iter<'this: 'node, 'node>(
        &'this self,
    ) -> impl use<'this, 'node, Node> + Iterator<Item = &'node Node>
    where
        Node: 'node,
    {
        Self::values(self)
    }
}

impl<Identifier, Node> Reset for FxHashMapSeed<Identifier, Node> {
    fn reset(&mut self) {
        self.clear();
    }
}
