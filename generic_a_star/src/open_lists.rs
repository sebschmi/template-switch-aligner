use binary_heap_plus::BinaryHeap;

use crate::{AStarNode, comparator::AStarNodeComparator, reset::Reset};

pub mod linear_heap;

/// The open list for the A* algorithm.
///
/// This is a priority queue that must sort nodes with [`AStarNodeComparator`].
pub trait AStarOpenList<Node: AStarNode>: Reset + Extend<Node> {
    /// Create a new empty open list.
    fn new() -> Self;

    /// Returns the number of nodes in the open list.
    fn len(&self) -> usize;

    /// Returns true if the open list is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Add the given node to the priority queue.
    fn push(&mut self, node: Node);

    /// Removes and returns a minimum element from the priority queue.
    ///
    /// Returns `None` if the priority queue is empty.
    fn pop_min(&mut self) -> Option<Node>;
}

impl<Node: AStarNode> AStarOpenList<Node> for BinaryHeap<Node, AStarNodeComparator> {
    fn new() -> Self {
        BinaryHeap::from_vec(Vec::new())
    }

    fn len(&self) -> usize {
        Self::len(self)
    }

    fn push(&mut self, node: Node) {
        Self::push(self, node);
    }

    fn pop_min(&mut self) -> Option<Node> {
        Self::pop(self)
    }
}

impl<Node, Comparator> Reset for BinaryHeap<Node, Comparator> {
    fn reset(&mut self) {
        self.clear();
    }
}
