use binary_heap_plus::BinaryHeap;

use crate::{AStarNode, AStarOpenList, comparator::AStarNodeComparator, reset::Reset};

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
