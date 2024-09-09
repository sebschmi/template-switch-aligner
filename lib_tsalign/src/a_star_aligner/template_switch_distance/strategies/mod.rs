use std::marker::PhantomData;

use node_ord::NodeOrdStrategy;

pub mod node_ord;

pub trait AlignmentStrategySelector: Eq + Clone + std::fmt::Debug {
    type NodeOrd: NodeOrdStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategies<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
}

pub trait AlignmentStrategy: Eq + Clone + std::fmt::Debug {
    fn create_root() -> Self;

    fn generate_successor(&self) -> Self;
}

impl<Selector: AlignmentStrategySelector> AlignmentStrategy for AlignmentStrategies<Selector> {
    fn create_root() -> Self {
        Self {
            node_ord_strategy: Selector::NodeOrd::create_root(),
        }
    }

    fn generate_successor(&self) -> Self {
        self.clone()
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategySelection<NodeOrd: NodeOrdStrategy> {
    phantom_data: PhantomData<NodeOrd>,
}

impl<NodeOrd: NodeOrdStrategy> AlignmentStrategySelector for AlignmentStrategySelection<NodeOrd> {
    type NodeOrd = NodeOrd;
}
