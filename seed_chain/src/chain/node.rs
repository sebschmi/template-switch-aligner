use generic_a_star::{AStarNode, cost::AStarCost};

use crate::seed::ChainingAnchor;

pub mod display;

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub enum Identifier {
    Root,
    Anchor { anchor: ChainingAnchor },
    Target,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Node<Cost> {
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

impl<Cost: AStarCost> AStarNode for Node<Cost> {
    type Identifier = Identifier;

    type EdgeType = EdgeType;

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

impl<Cost: AStarCost> Node<Cost> {
    pub fn new_root() -> Self {
        Self {
            identifier: Identifier::Root,
            predecessor: None,
            cost: Cost::zero(),
        }
    }

    pub fn generate_successor(
        &self,
        successor_identifier: Identifier,
        cost_increment: Cost,
    ) -> Option<Self> {
        if cost_increment == Cost::max_value() {
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
