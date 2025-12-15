use rustc_hash::{FxHashMapSeed, FxSeededState};

use crate::{AStarIdentifier, AStarNode, reset::Reset};

/// The closed list for the A* algorithm.
pub trait AStarClosedList<Identifier: AStarIdentifier, Node: AStarNode>: Reset {
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
    /// If there was a node previously mapped to this identifier, then it is returned.
    /// Otherwise, `None` is returned.
    fn insert(&mut self, identifier: Identifier, node: Node) -> Option<Node>;

    /// Return a reference to the node specified by `identifier`.
    ///
    /// Returns `None` if no node with the given identifier exists.
    fn get(&self, identifier: &Identifier) -> Option<&Node>;

    /// Returns true if the given identifier is mapped to a node.
    fn contains_identifier(&self, identifier: &Identifier) -> bool {
        self.get(identifier).is_some()
    }
}

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
