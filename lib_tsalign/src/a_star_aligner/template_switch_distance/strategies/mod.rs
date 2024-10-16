use std::marker::PhantomData;

use compact_genome::interface::alphabet::Alphabet;
use node_ord::NodeOrdStrategy;
use template_switch_min_length::TemplateSwitchMinLengthStrategy;

use super::Context;

pub mod node_ord;
pub mod template_switch_min_length;

pub trait AlignmentStrategySelector: Eq + Clone + std::fmt::Debug {
    type Alphabet: Alphabet;
    type NodeOrd: NodeOrdStrategy;
    type TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategies<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
    pub template_switch_min_length_strategy: Selector::TemplateSwitchMinLength,
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
            template_switch_min_length_strategy: Selector::TemplateSwitchMinLength::create_root(
                context,
            ),
        }
    }

    fn generate_successor<Strategies: AlignmentStrategySelector>(
        &self,
        context: &Context<Strategies>,
    ) -> Self {
        Self {
            node_ord_strategy: self.node_ord_strategy.generate_successor(context),
            template_switch_min_length_strategy: self
                .template_switch_min_length_strategy
                .generate_successor(context),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategySelection<
    AlphabetType: Alphabet,
    NodeOrd: NodeOrdStrategy,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
> {
    phantom_data: PhantomData<(AlphabetType, NodeOrd, TemplateSwitchMinLength)>,
}

impl<
        AlphabetType: Alphabet + std::fmt::Debug + Clone + Eq,
        NodeOrd: NodeOrdStrategy,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
    > AlignmentStrategySelector
    for AlignmentStrategySelection<AlphabetType, NodeOrd, TemplateSwitchMinLength>
{
    type Alphabet = AlphabetType;
    type NodeOrd = NodeOrd;
    type TemplateSwitchMinLength = TemplateSwitchMinLength;
}
