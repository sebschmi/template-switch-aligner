use super::strategies::{AlignmentStrategiesNodeIdentifier, AlignmentStrategySelector};
use std::hash::Hash;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd)]
pub struct Identifier<Strategies: AlignmentStrategySelector> {
    pub kind: IdentifierKind,
    pub strategies: AlignmentStrategiesNodeIdentifier<Strategies>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub enum IdentifierKind {
    Primary {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
        /// Positive for left flank, negative for right flank.
        flank_index: isize,
    },
    PrimaryReentry {
        reference_index: usize,
        query_index: usize,
        gap_type: GapType,
        /// Positive for left flank, negative for right flank.
        flank_index: isize,
    },
    TemplateSwitchEntrance {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        template_switch_first_offset: isize,
    },
    Secondary {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        length: usize,
        /// The index that does not jump.
        primary_index: usize,
        /// The index that jumps.
        secondary_index: usize,
        gap_type: GapType,
    },
    TemplateSwitchExit {
        entrance_reference_index: usize,
        entrance_query_index: usize,
        template_switch_primary: TemplateSwitchPrimary,
        template_switch_secondary: TemplateSwitchSecondary,
        /// The index that does not jump.
        primary_index: usize,
        length_difference: isize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub enum GapType {
    Insertion,
    Deletion,
    None,
}

/// The primary sequence is the sequence for which the template switch does not jump.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TemplateSwitchPrimary {
    Reference,
    Query,
}

/// The secondary sequence is the sequence for which the template switch jumps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TemplateSwitchSecondary {
    Reference,
    Query,
}

impl IdentifierKind {
    pub const fn new_primary(
        reference_index: usize,
        query_index: usize,
        flank_index: isize,
        gap_type: GapType,
    ) -> Self {
        debug_assert!(reference_index != usize::MAX);
        debug_assert!(query_index != usize::MAX);
        debug_assert!(reference_index < isize::MAX as usize);
        debug_assert!(query_index < isize::MAX as usize);

        Self::Primary {
            reference_index,
            query_index,
            flank_index,
            gap_type,
        }
    }

    pub fn generate_primary_diagonal_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index: reference_index + 1,
                    query_index: query_index + 1,
                    flank_index,
                    gap_type: GapType::None,
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_primary_deletion_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index: reference_index + 1,
                    query_index,
                    flank_index,
                    gap_type: GapType::Deletion,
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_primary_insertion_successor(self, flank_index: isize) -> Self {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                Self::Primary {
                    reference_index,
                    query_index: query_index + 1,
                    flank_index,
                    gap_type: GapType::Insertion,
                }
            }
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_initial_template_switch_entrance_successors(
        self,
    ) -> impl Iterator<Item = Self> {
        match self {
            Self::Primary {
                reference_index: entrance_reference_index,
                query_index: entrance_query_index,
                ..
            }
            | Self::PrimaryReentry {
                reference_index: entrance_reference_index,
                query_index: entrance_query_index,
                ..
            } => {
                let template_switch_first_offset = 0;

                [
                    TemplateSwitchPrimary::Reference,
                    TemplateSwitchPrimary::Reference,
                    TemplateSwitchPrimary::Query,
                    TemplateSwitchPrimary::Query,
                ]
                .into_iter()
                .zip([
                    TemplateSwitchSecondary::Reference,
                    TemplateSwitchSecondary::Query,
                    TemplateSwitchSecondary::Reference,
                    TemplateSwitchSecondary::Query,
                ])
                .map(
                    move |(template_switch_primary, template_switch_secondary)| {
                        Self::TemplateSwitchEntrance {
                            entrance_reference_index,
                            entrance_query_index,
                            template_switch_primary,
                            template_switch_secondary,
                            template_switch_first_offset,
                        }
                    },
                )
            }

            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    pub fn generate_secondary_diagonal_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: length + 1,
                primary_index: primary_index + 1,
                secondary_index: secondary_index - 1,
                gap_type: GapType::None,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The secondary contains a base missing in the primary.
    pub fn generate_secondary_deletion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index: secondary_index - 1,
                gap_type: GapType::Deletion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// The secondary misses a base present in the primary.
    pub fn generate_secondary_insertion_successor(self) -> Self {
        match self {
            Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length,
                primary_index,
                secondary_index,
                ..
            } => Self::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: length + 1,
                primary_index: primary_index + 1,
                secondary_index,
                gap_type: GapType::Insertion,
            },
            other => unreachable!(
                "Function is only called on primary identifiers, but this is: {other}."
            ),
        }
    }

    /// Returns the anti-diagonal for variants where it exists, or [`usize::MAX`](core::primitive::usize::MAX) otherwise.
    pub const fn anti_diagonal(self) -> usize {
        match self {
            Self::Primary {
                reference_index,
                query_index,
                ..
            } => {
                debug_assert!(reference_index != usize::MAX);
                debug_assert!(query_index != usize::MAX);
                debug_assert!(reference_index < isize::MAX as usize);
                debug_assert!(query_index < isize::MAX as usize);

                reference_index + query_index
            }
            _ => usize::MAX,
        }
    }
}

impl<Strategies: AlignmentStrategySelector> Hash for Identifier<Strategies> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.kind.hash(state);
        self.strategies.hash(state);
    }
}
