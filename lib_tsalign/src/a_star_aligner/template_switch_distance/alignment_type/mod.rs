use crate::a_star_aligner::alignment_result::IAlignmentType;

use super::identifier::{TemplateSwitchPrimary, TemplateSwitchSecondary};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AlignmentType {
    /// The query contains a base that is missing from the reference.
    Insertion,
    /// The query is missing a base present in the reference.
    Deletion,
    /// The query contains a different base than the reference.
    Substitution,
    /// The query contains the same base as the reference.
    Match,
    /// A template switch entrance.
    TemplateSwitchEntrance {
        primary: TemplateSwitchPrimary,
        secondary: TemplateSwitchSecondary,
        first_offset: isize,
    },
    /// A template switch exit.
    TemplateSwitchExit { length_difference: isize },
    /// This node is the root node, hence it was not generated via alignment.
    Root,
    /// The root node of a secondary graph.
    SecondaryRoot,
    /// A reentry node into the primary graph, treated like a root.
    PrimaryReentry,
    /// A shortcut in the primary matrix.
    ///
    /// Used only for computing lower bounds.
    PrimaryShortcut {
        delta_reference: isize,
        delta_query: isize,
    },
}

impl IAlignmentType for AlignmentType {
    fn is_repeatable(&self) -> bool {
        match self {
            Self::Insertion
            | Self::Deletion
            | Self::Substitution
            | Self::Match
            | Self::Root
            | Self::SecondaryRoot
            | Self::PrimaryReentry => true,
            Self::TemplateSwitchEntrance { .. }
            | Self::TemplateSwitchExit { .. }
            | Self::PrimaryShortcut { .. } => false,
        }
    }

    fn is_repeated(&self, previous: &Self) -> bool {
        match (self, previous) {
            (
                Self::TemplateSwitchEntrance {
                    primary: primary_a,
                    secondary: secondary_a,
                    ..
                },
                Self::TemplateSwitchEntrance {
                    primary: primary_b,
                    secondary: secondary_b,
                    ..
                },
            ) => primary_a == primary_b && secondary_a == secondary_b,
            (Self::TemplateSwitchExit { .. }, Self::TemplateSwitchExit { .. }) => true,
            (Self::PrimaryShortcut { .. }, Self::PrimaryShortcut { .. }) => false,
            (a, b) => a == b,
        }
    }

    fn is_internal(&self) -> bool {
        matches!(
            self,
            Self::Root | Self::SecondaryRoot | Self::PrimaryReentry
        )
    }
}
