use std::collections::HashMap;

use binary_heap_plus::{BinaryHeap, MinComparator};

use crate::a_star_aligner::template_switch_distance::Node;
use crate::a_star_aligner::AlignmentGraphNode;
use crate::config::TemplateSwitchConfig;

use super::strategies::AlignmentStrategySelector;

pub struct Context<Strategies: AlignmentStrategySelector> {
    pub config: TemplateSwitchConfig<Strategies::Alphabet>,

    /// Can be used by strategies.
    pub open_list: BinaryHeap<Node<Strategies>, MinComparator>,
    /// Can be used by strategies.
    pub closed_list: HashMap<
        <Node<Strategies> as AlignmentGraphNode<Strategies::Alphabet>>::Identifier,
        Node<Strategies>,
    >,
}

impl<Strategies: AlignmentStrategySelector> Context<Strategies> {
    pub fn new(config: TemplateSwitchConfig<Strategies::Alphabet>) -> Self {
        Self {
            config,
            open_list: BinaryHeap::new_min(),
            closed_list: Default::default(),
        }
    }
}
