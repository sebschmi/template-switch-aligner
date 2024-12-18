use std::fmt::Display;

use compact_genome::interface::sequence::GenomeSequence;
use generic_a_star::AStarNode;
use identifier::{GapType, IdentifierKind, TemplateSwitchPrimary, TemplateSwitchSecondary};
use num_traits::SaturatingSub;
use strategies::{
    node_ord::NodeOrdStrategy, template_switch_min_length::TemplateSwitchMinLengthStrategy,
    AlignmentStrategiesNodeIdentifier, AlignmentStrategiesNodeMemory, AlignmentStrategySelector,
};

use crate::costs::cost::Cost;

mod alignment_type;
pub mod context;
pub mod display;
mod identifier;
pub mod lower_bounds;
pub mod strategies;

pub use alignment_type::AlignmentType;
pub use context::Context;
pub use identifier::Identifier;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Node<Strategies: AlignmentStrategySelector> {
    identifier: Identifier<Strategies>,
    predecessor: Option<Identifier<Strategies>>,
    predecessor_edge_type: AlignmentType,
    cost: Cost,
    a_star_lower_bound: Cost,
    strategies: AlignmentStrategiesNodeMemory<Strategies>,
}

impl<Strategies: AlignmentStrategySelector> AStarNode for Node<Strategies> {
    type Identifier = Identifier<Strategies>;

    type EdgeType = AlignmentType;

    fn identifier(&self) -> &Self::Identifier {
        &self.identifier
    }

    fn cost(&self) -> Cost {
        self.cost
    }

    fn a_star_lower_bound(&self) -> Cost {
        self.a_star_lower_bound
    }

    fn predecessor(&self) -> Option<&Self::Identifier> {
        self.predecessor.as_ref()
    }

    fn predecessor_edge_type(&self) -> Option<Self::EdgeType> {
        Some(self.predecessor_edge_type)
    }
}

impl<Strategies: AlignmentStrategySelector> Node<Strategies> {
    pub fn new_root_at<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        reference_index: usize,
        query_index: usize,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            identifier: Identifier {
                kind: IdentifierKind::Primary {
                    reference_index,
                    query_index,
                    gap_type: GapType::None,
                    flank_index: 0,
                },
                strategies: AlignmentStrategiesNodeIdentifier::create_root(context),
            },
            predecessor: None,
            predecessor_edge_type: AlignmentType::Root,
            cost: Cost::ZERO,
            a_star_lower_bound: Cost::ZERO,
            strategies: AlignmentStrategiesNodeMemory::create_root(context),
        }
    }

    fn generate_primary_diagonal_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        is_match: bool,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ (IdentifierKind::Primary { flank_index, .. }
        | IdentifierKind::PrimaryReentry { flank_index, .. }) = self.identifier.kind
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_diagonal_successor(successor_flank_index),
            cost_increment,
            match (
                is_match,
                flank_index == successor_flank_index && successor_flank_index == 0,
            ) {
                (true, true) => AlignmentType::PrimaryMatch,
                (true, false) => AlignmentType::PrimaryFlankMatch,
                (false, true) => AlignmentType::PrimarySubstitution,
                (false, false) => AlignmentType::PrimaryFlankSubstitution,
            },
            context,
        ))
    }

    fn generate_primary_deletion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier {
            kind:
                IdentifierKind::Primary { flank_index, .. }
                | IdentifierKind::PrimaryReentry { flank_index, .. },
            ..
        } = self.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_deletion_successor(successor_flank_index),
            cost_increment,
            if flank_index == successor_flank_index && successor_flank_index == 0 {
                AlignmentType::PrimaryDeletion
            } else {
                AlignmentType::PrimaryFlankDeletion
            },
            context,
        ))
    }

    fn generate_primary_insertion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        successor_flank_index: isize,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ (IdentifierKind::Primary { flank_index, .. }
        | IdentifierKind::PrimaryReentry { flank_index, .. }) = self.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_primary_insertion_successor(successor_flank_index),
            cost_increment,
            if flank_index == successor_flank_index && successor_flank_index == 0 {
                AlignmentType::PrimaryInsertion
            } else {
                AlignmentType::PrimaryFlankInsertion
            },
            context,
        ))
    }

    fn generate_initial_template_switch_entrance_successors<
        'result,
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &'result self,
        cost_increment: Cost,
        context: &'result Context<SubsequenceType, Strategies>,
    ) -> impl 'result + Iterator<Item = Self> {
        if !matches!(
            self.node_data.identifier,
            Identifier::Primary { .. } | Identifier::PrimaryReentry { .. }
        ) {
            unreachable!("This method is only called on primary nodes.")
        }

        self.node_data
            .identifier
            .generate_initial_template_switch_entrance_successors()
            .map(move |identifier| {
                let Identifier::TemplateSwitchEntrance {
                    template_switch_primary,
                    template_switch_secondary,
                    template_switch_first_offset,
                    ..
                } = &identifier
                else {
                    unreachable!("This closure is only called on template switch entrances.")
                };

                self.generate_successor(
                    identifier,
                    cost_increment,
                    AlignmentType::TemplateSwitchEntrance {
                        primary: *template_switch_primary,
                        secondary: *template_switch_secondary,
                        first_offset: *template_switch_first_offset,
                    },
                    context,
                )
            })
    }

    fn generate_template_switch_entrance_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        successor_template_switch_first_offset: isize,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchEntrance {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                template_switch_first_offset: successor_template_switch_first_offset,
            },
            cost_increment,
            AlignmentType::TemplateSwitchEntrance {
                primary: template_switch_primary,
                secondary: template_switch_secondary,
                first_offset: successor_template_switch_first_offset,
            },
            context,
        ))
    }

    fn generate_secondary_root_node<
        'this,
        'reference,
        'query,
        'context,
        SubsequenceType: compact_genome::interface::sequence::GenomeSequence<
                Strategies::Alphabet,
                SubsequenceType,
            > + ?Sized,
    >(
        &'this self,
        context: &'context mut Context<'reference, 'query, SubsequenceType, Strategies>,
    ) -> impl use<'this, 'reference, 'query, 'context, SubsequenceType, Strategies>
           + IntoIterator<Item = Self> {
        let Identifier::TemplateSwitchEntrance {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            template_switch_first_offset,
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch entrance nodes.")
        };

        let primary_index = match template_switch_primary {
            TemplateSwitchPrimary::Reference => entrance_reference_index,
            TemplateSwitchPrimary::Query => entrance_query_index,
        };

        let secondary_index = (match template_switch_secondary {
            TemplateSwitchSecondary::Reference => entrance_reference_index,
            TemplateSwitchSecondary::Query => entrance_query_index,
        } as isize
            + template_switch_first_offset) as usize;

        let secondary_root_node = self.generate_successor(
            Identifier::Secondary {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                length: 0,
                primary_index,
                secondary_index,
                gap_type: GapType::None,
            },
            0.into(),
            AlignmentType::SecondaryRoot,
            context,
        );

        self.strategies
            .template_switch_min_length_strategy
            .template_switch_min_length_lookahead(secondary_root_node, context)
    }

    fn generate_secondary_diagonal_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        is_match: bool,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_diagonal_successor(),
            cost_increment,
            if is_match {
                AlignmentType::SecondaryMatch
            } else {
                AlignmentType::SecondarySubstitution
            },
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_deletion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_deletion_successor(),
            cost_increment,
            AlignmentType::SecondaryDeletion,
            context,
        ))
    }

    /// The secondary contains a base missing in the primary.
    fn generate_secondary_insertion_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let predecessor_identifier @ Identifier::Secondary { .. } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            predecessor_identifier.generate_secondary_insertion_successor(),
            cost_increment,
            AlignmentType::SecondaryInsertion,
            context,
        ))
    }

    fn generate_initial_template_switch_exit_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::Secondary {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            primary_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on secondary nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                length_difference: 0,
            },
            cost_increment,
            AlignmentType::TemplateSwitchExit {
                length_difference: 0,
            },
            context,
        ))
    }

    fn generate_template_switch_exit_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        cost_increment: Cost,
        successor_length_difference: isize,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            template_switch_secondary,
            primary_index,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        Some(self.generate_successor(
            Identifier::TemplateSwitchExit {
                entrance_reference_index,
                entrance_query_index,
                template_switch_primary,
                template_switch_secondary,
                primary_index,
                length_difference: successor_length_difference,
            },
            cost_increment,
            AlignmentType::TemplateSwitchExit {
                length_difference: successor_length_difference,
            },
            context,
        ))
    }

    pub fn generate_primary_reentry_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        let Identifier::TemplateSwitchExit {
            entrance_reference_index,
            entrance_query_index,
            template_switch_primary,
            primary_index,
            length_difference,
            ..
        } = self.node_data.identifier
        else {
            unreachable!("This method is only called on template switch exit nodes.")
        };

        let (reference_index, query_index) = match template_switch_primary {
            TemplateSwitchPrimary::Reference => {
                let primary_length = primary_index - entrance_reference_index;
                let anti_primary_length = primary_length as isize + length_difference;
                let query_index = entrance_query_index as isize + anti_primary_length;

                if query_index < 0 {
                    return None;
                }

                (primary_index, query_index as usize)
            }
            TemplateSwitchPrimary::Query => {
                let primary_length = primary_index - entrance_query_index;
                let anti_primary_length = primary_length as isize + length_difference;
                let reference_index = entrance_reference_index as isize + anti_primary_length;

                if reference_index < 0 {
                    return None;
                }

                (reference_index as usize, primary_index)
            }
        };

        debug_assert!(reference_index != usize::MAX, "{self:?}");
        debug_assert!(query_index != usize::MAX, "{self:?}");
        debug_assert!(reference_index < isize::MAX as usize, "{self:?}");
        debug_assert!(query_index < isize::MAX as usize, "{self:?}");

        Some(self.generate_successor(
            Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: -context.config.right_flank_length,
            },
            0.into(),
            AlignmentType::PrimaryReentry,
            context,
        ))
    }

    fn generate_primary_shortcut_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        delta_reference: isize,
        delta_query: isize,
        cost_increment: Cost,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Option<Self> {
        if cost_increment == Cost::MAX {
            return None;
        }

        let (Identifier::Primary {
            reference_index,
            query_index,
            flank_index,
            ..
        }
        | Identifier::PrimaryReentry {
            reference_index,
            query_index,
            flank_index,
            ..
        }) = self.node_data.identifier
        else {
            unreachable!("This method is only called on primary nodes.")
        };
        assert_eq!(flank_index, context.config.left_flank_length);

        let reference_index =
            usize::try_from(isize::try_from(reference_index).unwrap() + delta_reference).ok()?;
        let query_index =
            usize::try_from(isize::try_from(query_index).unwrap() + delta_query).ok()?;

        if reference_index >= context.reference.len() || query_index >= context.query.len() {
            return None;
        }

        Some(self.generate_successor(
            Identifier::PrimaryReentry {
                reference_index,
                query_index,
                gap_type: GapType::None,
                flank_index: -context.config.right_flank_length,
            },
            cost_increment,
            AlignmentType::PrimaryShortcut {
                delta_reference,
                delta_query,
            },
            context,
        ))
    }

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        identifier: Identifier<Strategies>,
        cost_increment: Cost,
        alignment_type: AlignmentType,
        context: &Context<SubsequenceType, Strategies>,
    ) -> Self {
        let cost = self.cost + cost_increment;
        let a_star_lower_bound = self.a_star_lower_bound.saturating_sub(&cost_increment);

        Self {
            identifier,
            predecessor: Some(self.identifier),
            predecessor_edge_type: alignment_type,
            cost,
            a_star_lower_bound,
            strategies: self
                .strategies
                .generate_successor(identifier, alignment_type, context),
        }
    }
}

impl NodeData {
    fn lower_bound_cost(&self) -> Cost {
        self.cost + self.a_star_lower_bound
    }
}

impl<Strategies: AlignmentStrategySelector> Ord for Node<Strategies> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.strategies.node_ord_strategy.cmp(self, other)
    }
}

impl<Strategies: AlignmentStrategySelector> PartialOrd for Node<Strategies> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Strategies: AlignmentStrategySelector> Display for Node<Strategies>
where
    AlignmentStrategiesNodeMemory<Strategies>: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            identifier,
            predecessor,
            predecessor_edge_type,
            cost,
            a_star_lower_bound,
            strategies,
        } = self;
        write!(f, "{identifier}; ")?;
        if let Some(predecessor) = predecessor {
            write!(f, "predecessor: {predecessor}; ")?;
        }
        write!(f, "alignment_type: {predecessor_edge_type}; ")?;
        write!(f, "cost: {cost} + {a_star_lower_bound}; ")?;
        write!(f, "strategies: {strategies}")
    }
}
