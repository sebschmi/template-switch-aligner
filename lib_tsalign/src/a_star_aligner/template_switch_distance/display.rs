use std::fmt::{Display, Formatter, Result};

use super::{AlignmentType, GapType, Identifier, TemplateSwitchPrimary, TemplateSwitchSecondary};

impl Display for AlignmentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::PrimaryInsertion | Self::PrimaryFlankInsertion | Self::SecondaryInsertion => {
                write!(f, "I")
            }
            Self::PrimaryDeletion | Self::PrimaryFlankDeletion | Self::SecondaryDeletion => {
                write!(f, "D")
            }
            Self::PrimarySubstitution
            | Self::PrimaryFlankSubstitution
            | Self::SecondarySubstitution => write!(f, "S"),
            Self::PrimaryMatch | Self::PrimaryFlankMatch | Self::SecondaryMatch => write!(f, "M"),
            Self::TemplateSwitchEntrance {
                primary,
                secondary,
                first_offset,
            } => write!(f, "[TS{primary}{secondary}{first_offset}:"),
            Self::TemplateSwitchExit { anti_primary_gap } => write!(f, ":{anti_primary_gap}]"),
            Self::Root => Ok(()),
            Self::SecondaryRoot => Ok(()),
            Self::PrimaryReentry => Ok(()),
            Self::PrimaryShortcut {
                delta_reference,
                delta_query,
            } => write!(f, "[PS:R{delta_reference}Q{delta_query}]"),
        }
    }
}

impl Display for GapType {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Insertion => write!(f, "I"),
            Self::Deletion => write!(f, "D"),
            Self::None => write!(f, "M/S"),
        }
    }
}

impl Display for TemplateSwitchPrimary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Reference => write!(f, "R"),
            Self::Query => write!(f, "Q"),
        }
    }
}

impl Display for TemplateSwitchSecondary {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Reference => write!(f, "R"),
            Self::Query => write!(f, "Q"),
        }
    }
}

impl<PrimaryExtraData> Display for Identifier<PrimaryExtraData> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                flank_index,
                gap_type,
                ..
            } => write!(
                f,
                "Primary({}R, {}Q, {}F, {})",
                reference_index, query_index, flank_index, gap_type
            ),

            Self::PrimaryReentry {
                reference_index,
                query_index,
                flank_index,
                gap_type,
                ..
            } => write!(
                f,
                "PrimaryReentry({}R, {}Q, {}F, {})",
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
                    "TemplateSwitchEntrance({}R, {}Q, {}P, {}S, {}O)",
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
                length,
                primary_index,
                secondary_index,
                gap_type,
            } => write!(
                f,
                "Secondary({}R, {}Q, {}L, {}P, {}S, {}, {}, {})",
                entrance_reference_index,
                entrance_query_index,
                length,
                primary_index,
                secondary_index,
                template_switch_primary,
                template_switch_secondary,
                gap_type
            ),

            Self::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                anti_primary_gap,
            } => write!(
                f,
                "TemplateSwitchExit({}R, {}Q, {}P, {}G, {}, {})",
                entrance_reference_index,
                entrance_query_index,
                primary_index,
                anti_primary_gap,
                template_switch_primary,
                template_switch_secondary
            ),
        }
    }
}
