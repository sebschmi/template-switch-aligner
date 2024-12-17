use generic_a_star::{cost::Cost, AStarNode};

use crate::seed::ChainingAnchor;

pub mod display;

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub enum Identifier {
    Root,
    Anchor { anchor: ChainingAnchor },
    Target,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Node {
    identifier: Identifier,
    predecessor: Option<Identifier>,
    cost: Cost,
}

#[derive(Debug, Clone)]
pub struct EdgeType {
    #[expect(dead_code)]
    pub from: Identifier,
    pub to: Identifier,
}

impl AStarNode for Node {
    type Identifier = Identifier;

    type EdgeType = EdgeType;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Cost {
        Cost::ZERO
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        self.predecessor.clone().map(|from| EdgeType {
            from,
            to: self.identifier.clone(),
        })
    }
}

impl Node {
    pub fn new_root() -> Self {
        Self {
            identifier: Identifier::Root,
            predecessor: None,
            cost: Cost::ZERO,
        }
    }

    pub fn generate_successor(
        &self,
        successor_identifier: Identifier,
        cost_increment: Cost,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        debug_assert!(!matches!(self.identifier, Identifier::Target));

        Some(Self {
            identifier: successor_identifier,
            predecessor: Some(self.identifier.clone()),
            cost: self.cost + cost_increment,
        })
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.cmp(&other.cost)
    }
}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
