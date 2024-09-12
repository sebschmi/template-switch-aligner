use std::fmt::{Display, Formatter, Result};

use super::{AlignmentType, GapType, Identifier, TemplateSwitchOrigin, TemplateSwitchTarget};

impl Display for AlignmentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            AlignmentType::Insertion => write!(f, "I"),
            AlignmentType::Deletion => write!(f, "D"),
            AlignmentType::Substitution => write!(f, "S"),
            AlignmentType::Match => write!(f, "M"),
            AlignmentType::TemplateSwitchEntrance {
                origin,
                target,
                first_offset,
            } => write!(f, "[TS{origin}{target}{first_offset}]"),
            AlignmentType::Root => Ok(()),
        }
    }
}

impl Display for GapType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            GapType::Insertion => write!(f, "I"),
            GapType::Deletion => write!(f, "D"),
            GapType::None => write!(f, "M/S"),
        }
    }
}

impl Display for TemplateSwitchOrigin {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            TemplateSwitchOrigin::Reference => write!(f, "R"),
            TemplateSwitchOrigin::Query => write!(f, "Q"),
        }
    }
}

impl Display for TemplateSwitchTarget {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            TemplateSwitchTarget::Reference => write!(f, "R"),
            TemplateSwitchTarget::Query => write!(f, "Q"),
        }
    }
}

impl Display for Identifier {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
            } => write!(
                f,
                "Primary({}R, {}Q, {}F, {})",
                reference_index, query_index, flank_index, gap_type
            ),

            Self::TemplateSwitchEntrance {
                entrance_reference_index: origin_reference_index,
                entrance_query_index: origin_query_index,
                template_switch_origin,
                template_switch_target,
                template_switch_first_offset,
            } => {
                write!(
                    f,
                    "TemplateSwitchEntrance({}R, {}Q, {}O, {}T, {}O1)",
                    origin_reference_index,
                    origin_query_index,
                    template_switch_origin,
                    template_switch_target,
                    template_switch_first_offset
                )
            }
        }
    }
}
