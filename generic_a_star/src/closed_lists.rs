use rustc_hash::{FxHashMapSeed, FxSeededState};

use crate::{AStarClosedList, AStarIdentifier, AStarNode, reset::Reset};

impl<Identifier: AStarIdentifier, Node: AStarNode> AStarClosedList<Identifier, Node>
    for FxHashMapSeed<Identifier, Node>
{
    fn new() -> Self {
        FxHashMapSeed::with_hasher(FxSeededState::with_seed(0))
    }

    fn len(&self) -> usize {
        Self::len(self)
    }

    fn insert(&mut self, identifier: Identifier, node: Node) -> Option<Node> {
        Self::insert(self, identifier, node)
    }

    fn get(&self, identifier: &Identifier) -> Option<&Node> {
        Self::get(self, identifier)
    }
}

impl<Identifier, Node> Reset for FxHashMapSeed<Identifier, Node> {
    fn reset(&mut self) {
        self.clear();
    }
}
