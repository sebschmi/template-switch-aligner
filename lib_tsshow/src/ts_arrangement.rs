use lib_tsalign::a_star_aligner::template_switch_distance::{
    AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary,
};
use strong_type::StrongType;

pub mod complement;
pub mod inner;
pub mod source;

pub struct TemplateSwitch {
    primary: TemplateSwitchPrimary,
    secondary: TemplateSwitchSecondary,
    sp1_reference: ArrangementCharColumn,
    sp1_query: ArrangementCharColumn,
    sp4_reference: ArrangementCharColumn,
    sp4_query: ArrangementCharColumn,
    sp2_secondary: SourceColumn,
    sp3_secondary: SourceColumn,
    inner_alignment: Vec<AlignmentType>,
}

/// A column index in a source string.
#[derive(StrongType)]
#[strong_type(conversion)]
pub struct SourceColumn(usize);

/// A column index over the characters in an arrangement.
#[derive(StrongType)]
#[strong_type(conversion)]
pub struct ArrangementCharColumn(usize);

/// A column index in an arrangement.
///
/// This is the index that represents where the character will actually be rendered.
#[derive(StrongType)]
#[strong_type(conversion)]
pub struct ArrangementColumn(usize);
