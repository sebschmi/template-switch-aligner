use std::marker::PhantomData;

use compact_genome::interface::alphabet::Alphabet;
use node_ord::NodeOrdStrategy;

use super::Context;

pub mod node_ord;
pub mod template_switch_min_length;

pub trait AlignmentStrategySelector: Eq + Clone + std::fmt::Debug {
    type Alphabet: Alphabet;
    type NodeOrd: NodeOrdStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategies<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
}

pub trait AlignmentStrategy: Eq + Clone + std::fmt::Debug {
    fn create_root<Strategies: AlignmentStrategySelector>(context: &Context<Strategies>) -> Self;

    fn generate_successor<Strategies: AlignmentStrategySelector>(
        &self,
        context: &Context<Strategies>,
    ) -> Self;
}

impl<Selector: AlignmentStrategySelector> AlignmentStrategy for AlignmentStrategies<Selector> {
    fn create_root<Strategies: AlignmentStrategySelector>(context: &Context<Strategies>) -> Self {
        Self {
            node_ord_strategy: Selector::NodeOrd::create_root(context),
        }
    }

    fn generate_successor<Strategies: AlignmentStrategySelector>(
        &self,
        context: &Context<Strategies>,
    ) -> Self {
        Self {
            node_ord_strategy: self.node_ord_strategy.generate_successor(context),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategySelection<AlphabetType: Alphabet, NodeOrd: NodeOrdStrategy> {
    phantom_data: PhantomData<(AlphabetType, NodeOrd)>,
}

impl<AlphabetType: Alphabet + std::fmt::Debug + Clone + Eq, NodeOrd: NodeOrdStrategy>
    AlignmentStrategySelector for AlignmentStrategySelection<AlphabetType, NodeOrd>
{
    type Alphabet = AlphabetType;
    type NodeOrd = NodeOrd;
}
