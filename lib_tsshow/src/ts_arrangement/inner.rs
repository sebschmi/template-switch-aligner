use super::{SourceColumn, TemplateSwitch, source::TsSourceArrangement};

pub struct TsInnerArrangement {
    reference_inners: Vec<Vec<InnerChar>>,
    query_inners: Vec<Vec<InnerChar>>,
}

#[derive(Debug, Clone, Copy)]
pub enum InnerChar {
    Inner {
        column: SourceColumn,
        lower_case: bool,
    },
    Gap,
    Blank,
}

impl TsInnerArrangement {
    pub fn new(
        source_arrangement: &TsSourceArrangement,
        template_switches: &[TemplateSwitch],
    ) -> Self {
        todo!()
    }
}
