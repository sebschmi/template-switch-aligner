use crate::config::TemplateSwitchConfig;

use super::strategies::AlignmentStrategySelector;

pub struct Context<Strategies: AlignmentStrategySelector> {
    pub config: TemplateSwitchConfig<Strategies::Alphabet>,
}

impl<Strategies: AlignmentStrategySelector> Context<Strategies> {
    pub fn new(config: TemplateSwitchConfig<Strategies::Alphabet>) -> Self {
        Self { config }
    }
}
