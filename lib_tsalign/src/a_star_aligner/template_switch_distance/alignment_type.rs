use crate::a_star_aligner::alignment_result::IAlignmentType;

use super::identifier::{TemplateSwitchPrimary, TemplateSwitchSecondary};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AlignmentType {
    /// The query contains a base that is missing from the reference.
    PrimaryInsertion,
    /// The query is missing a base present in the reference.
    PrimaryDeletion,
    /// The query contains a different base than the reference.
    PrimarySubstitution,
    /// The query contains the same base as the reference.
    PrimaryMatch,
    /// The query contains a base that is missing from the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankInsertion,
    /// The query is missing a base present in the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankDeletion,
    /// The query contains a different base than the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankSubstitution,
    /// The query contains the same base as the reference.
    ///
    /// This happens inside a TS flank.
    PrimaryFlankMatch,
    /// The TS secondary contains a base that is missing from the TS primary.
    SecondaryInsertion,
    /// The TS secondary is missing a base present in the TS primary.
    SecondaryDeletion,
    /// The TS secondary contains a different base than the TS primary.
    SecondarySubstitution,
    /// The TS secondary contains the same base as the TS primary.
    SecondaryMatch,
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
            Self::PrimaryInsertion
            | Self::PrimaryFlankInsertion
            | Self::SecondaryInsertion
            | Self::PrimaryDeletion
            | Self::PrimaryFlankDeletion
            | Self::SecondaryDeletion
            | Self::PrimarySubstitution
            | Self::PrimaryFlankSubstitution
            | Self::SecondarySubstitution
            | Self::PrimaryMatch
            | Self::PrimaryFlankMatch
            | Self::SecondaryMatch
            | Self::Root
            | Self::SecondaryRoot
            | Self::PrimaryReentry => true,
            Self::TemplateSwitchEntrance { .. }
            | Self::TemplateSwitchExit { .. }
            | Self::PrimaryShortcut { .. } => false,
        }
    }

    fn is_repeated(&self, previous: &Self) -> bool {
        if (matches!(self, Self::PrimaryInsertion | Self::PrimaryFlankInsertion)
            && matches!(
                previous,
                Self::PrimaryInsertion | Self::PrimaryFlankInsertion
            ))
            || (matches!(self, Self::PrimaryDeletion | Self::PrimaryFlankDeletion)
                && matches!(previous, Self::PrimaryDeletion | Self::PrimaryFlankDeletion))
            || (matches!(
                self,
                Self::PrimarySubstitution | Self::PrimaryFlankSubstitution
            ) && matches!(
                previous,
                Self::PrimarySubstitution | Self::PrimaryFlankSubstitution
            ))
            || (matches!(self, Self::PrimaryMatch | Self::PrimaryFlankMatch)
                && matches!(previous, Self::PrimaryMatch | Self::PrimaryFlankMatch))
        {
            return true;
        }

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
