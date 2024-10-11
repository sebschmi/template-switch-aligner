use crate::config::TemplateSwitchConfig;

pub struct Context<AlphabetType> {
    pub config: TemplateSwitchConfig<AlphabetType>,
}

impl<AlphabetType> Context<AlphabetType> {
    pub fn new(config: TemplateSwitchConfig<AlphabetType>) -> Self {
        Self { config }
    }
}
