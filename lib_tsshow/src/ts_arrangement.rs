use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;
use strong_type::StrongType;
use tagged_vec::TaggedVec;

pub mod source;

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
