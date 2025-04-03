use super::{
    index_types::{ArrangementCharColumn, SourceColumn},
    source::SourceChar,
};
use lib_tsalign::a_star_aligner::{
    alignment_result::alignment::Alignment,
    template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
};

#[derive(Debug, Clone)]
pub struct TemplateSwitch {
    pub primary: TemplateSwitchPrimary,
    pub secondary: TemplateSwitchSecondary,
    pub sp1_reference: ArrangementCharColumn,
    pub sp1_query: ArrangementCharColumn,
    pub sp4_reference: ArrangementCharColumn,
    pub sp4_query: ArrangementCharColumn,
    pub sp2_secondary: SourceColumn,
    pub sp3_secondary: SourceColumn,
    pub inner: Vec<SourceChar>,
    pub inner_alignment: Alignment<AlignmentType>,
}
