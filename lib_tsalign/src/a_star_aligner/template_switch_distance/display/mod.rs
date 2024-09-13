use std::fmt::{Display, Formatter, Result};

use super::{AlignmentType, GapType, Identifier, TemplateSwitchPrimary, TemplateSwitchSecondary};

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

impl Display for TemplateSwitchPrimary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            TemplateSwitchPrimary::Reference => write!(f, "R"),
            TemplateSwitchPrimary::Query => write!(f, "Q"),
        }
    }
}

impl Display for TemplateSwitchSecondary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            TemplateSwitchSecondary::Reference => write!(f, "R"),
            TemplateSwitchSecondary::Query => write!(f, "Q"),
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
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_first_offset,
            } => {
                write!(
                    f,
                    "TemplateSwitchEntrance({}R, {}Q, {}P, {}S, {}O1)",
                    entrance_reference_index,
                    entrance_query_index,
                    template_switch_primary,
                    template_switch_secondary,
                    template_switch_first_offset
                )
            }

            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                secondary_index,
                gap_type,
            } => write!(
                f,
                "Secondary({}R, {}Q, {}P, {}S, {}, {}, {})",
                entrance_reference_index,
                entrance_query_index,
                primary_index,
                secondary_index,
                template_switch_primary,
                template_switch_secondary,
                gap_type
            ),
        }
    }
}
