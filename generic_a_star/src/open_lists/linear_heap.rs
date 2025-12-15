use std::collections::VecDeque;

use num_traits::Zero;

use crate::{AStarNode, cost::AStarCost, open_lists::AStarOpenList, reset::Reset};

pub struct LinearHeap<Node: AStarNode> {
    heap: VecDeque<Vec<Node>>,
    cost_offset: Node::Cost,
    len: usize,
}

impl<Node: AStarNode> AStarOpenList<Node> for LinearHeap<Node> {
    fn new() -> Self {
        Self {
            heap: Default::default(),
            cost_offset: Zero::zero(),
            len: 0,
        }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn push(&mut self, node: Node) {
        let cost = node.cost();

        if cost < self.cost_offset {
            let prefix_len = (self.cost_offset - cost).as_usize();
            self.heap.reserve(prefix_len);
            for _ in 0..prefix_len {
                self.heap.push_front(Vec::new());
            }
            self.cost_offset = cost;

            self.heap[0].push(node);
        } else if cost >= self.cost_offset + Node::Cost::from_usize(self.heap.len()) {
            let suffix_len = (cost + Node::Cost::from(1u8)
                - (self.cost_offset + Node::Cost::from_usize(self.heap.len())))
            .as_usize();
            self.heap.reserve(suffix_len);
            for _ in 0..suffix_len {
                self.heap.push_back(Vec::new());
            }

            self.heap.back_mut().unwrap().push(node);
        } else {
            let index = (cost - self.cost_offset).as_usize();
            self.heap[index].push(node);
        }

        self.len += 1;
    }

    fn pop_min(&mut self) -> Option<Node> {
        while let Some(front) = self.heap.front_mut() {
            if front.is_empty() {
                // Lazily pop the front of the heap.
                self.heap.pop_front();
                self.cost_offset += Node::Cost::from(1u8);
            } else {
                self.len -= 1;
                return front.pop();
            }
        }

        None
    }
}

impl<Node: AStarNode> Extend<Node> for LinearHeap<Node> {
    fn extend<T: IntoIterator<Item = Node>>(&mut self, iter: T) {
        for node in iter {
            self.push(node);
        }
    }
}

impl<Node: AStarNode> Reset for LinearHeap<Node> {
    fn reset(&mut self) {
        self.heap.clear();
        self.cost_offset = Zero::zero();
        self.len = 0;
    }
}
