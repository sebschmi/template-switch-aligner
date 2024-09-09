use std::marker::PhantomData;

use node_ord::NodeOrdStrategy;

use super::Context;

pub mod flank_cost;
pub mod node_ord;

pub trait AlignmentStrategySelector: Eq + Clone + std::fmt::Debug {
    type NodeOrd: NodeOrdStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategies<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
}

pub trait AlignmentStrategy: Eq + Clone + std::fmt::Debug {
    fn create_root(context: &Context) -> Self;

    fn generate_successor(&self, context: &Context) -> Self;
}

impl<Selector: AlignmentStrategySelector> AlignmentStrategy for AlignmentStrategies<Selector> {
    fn create_root(context: &Context) -> Self {
        Self {
            node_ord_strategy: Selector::NodeOrd::create_root(context),
        }
    }

    fn generate_successor(&self, context: &Context) -> Self {
        Self {
            node_ord_strategy: self.node_ord_strategy.generate_successor(context),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategySelection<NodeOrd: NodeOrdStrategy> {
    phantom_data: PhantomData<NodeOrd>,
}

impl<NodeOrd: NodeOrdStrategy> AlignmentStrategySelector for AlignmentStrategySelection<NodeOrd> {
    type NodeOrd = NodeOrd;
}
