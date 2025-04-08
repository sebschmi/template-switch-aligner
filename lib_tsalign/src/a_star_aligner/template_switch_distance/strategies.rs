use std::{fmt::Debug, marker::PhantomData};

use chaining::ChainingStrategy;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::cost::AStarCost;
use node_ord::NodeOrdStrategy;
use primary_match::PrimaryMatchStrategy;
use primary_range::PrimaryRangeStrategy;
use secondary_deletion::SecondaryDeletionStrategy;
use shortcut::ShortcutStrategy;
use template_switch_count::TemplateSwitchCountStrategy;
use template_switch_min_length::TemplateSwitchMinLengthStrategy;

use super::{AlignmentType, Context, Identifier};

pub mod chaining;
pub mod node_ord;
pub mod primary_match;
pub mod primary_range;
pub mod secondary_deletion;
pub mod shortcut;
pub mod template_switch_count;
pub mod template_switch_min_length;

pub trait AlignmentStrategySelector: Eq + Clone + std::fmt::Debug {
    type Alphabet: Alphabet;
    type Cost: AStarCost;
    type NodeOrd: NodeOrdStrategy<Self::Cost, Self::PrimaryMatch>;
    type TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Self::Cost>;
    type Chaining: ChainingStrategy<Self::Cost>;
    type TemplateSwitchCount: TemplateSwitchCountStrategy;
    type SecondaryDeletion: SecondaryDeletionStrategy;
    type Shortcut: ShortcutStrategy<Self::Cost>;
    type PrimaryMatch: PrimaryMatchStrategy<Self::Cost>;
    type PrimaryRange: PrimaryRangeStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategiesNodeMemory<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
    pub template_switch_min_length_strategy: Selector::TemplateSwitchMinLength,
    pub template_switch_count: Selector::TemplateSwitchCount,
    pub primary_match: Selector::PrimaryMatch,
}

pub trait AlignmentStrategy: Eq + Clone + std::fmt::Debug {
    fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;

    fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
        Strategies: AlignmentStrategySelector,
    >(
        &self,
        identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self;
}

impl<Strategies: AlignmentStrategySelector> AlignmentStrategiesNodeMemory<Strategies> {
    pub fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            node_ord_strategy: Strategies::NodeOrd::create_root(context),
            template_switch_min_length_strategy: Strategies::TemplateSwitchMinLength::create_root(
                context,
            ),
            template_switch_count: Strategies::TemplateSwitchCount::create_root(context),
            primary_match: Strategies::PrimaryMatch::create_root(context),
        }
    }

    pub fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        identifier: Identifier<
            <<Strategies as AlignmentStrategySelector>::PrimaryMatch as PrimaryMatchStrategy<
                <Strategies as AlignmentStrategySelector>::Cost,
            >>::IdentifierPrimaryExtraData,
        >,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            node_ord_strategy: self.node_ord_strategy.generate_successor(
                identifier,
                alignment_type,
                context,
            ),
            template_switch_min_length_strategy: self
                .template_switch_min_length_strategy
                .generate_successor(identifier, alignment_type, context),
            template_switch_count: self.template_switch_count.generate_successor(
                identifier,
                alignment_type,
                context,
            ),
            primary_match: self.primary_match.generate_successor(
                identifier,
                alignment_type,
                context,
            ),
        }
    }
}

pub struct AlignmentStrategySelection<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> {
    #[allow(clippy::type_complexity)]
    phantom_data: PhantomData<(
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    )>,
}

impl<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> AlignmentStrategySelector
    for AlignmentStrategySelection<
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    >
{
    type Alphabet = AlphabetType;
    type Cost = Cost;
    type NodeOrd = NodeOrd;
    type TemplateSwitchMinLength = TemplateSwitchMinLength;
    type Chaining = Chaining;
    type TemplateSwitchCount = TemplateSwitchCount;
    type SecondaryDeletion = SecondaryDeletion;
    type Shortcut = Shortcut;
    type PrimaryMatch = PrimaryMatch;
    type PrimaryRange = PrimaryRange;
}

impl<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> Debug
    for AlignmentStrategySelection<
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    >
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentStrategySelection").finish()
    }
}

impl<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> Clone
    for AlignmentStrategySelection<
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    >
{
    fn clone(&self) -> Self {
        Self {
            phantom_data: self.phantom_data,
        }
    }
}

impl<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> PartialEq
    for AlignmentStrategySelection<
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    >
{
    fn eq(&self, other: &Self) -> bool {
        self.phantom_data == other.phantom_data
    }
}

impl<
    AlphabetType: Alphabet,
    Cost: AStarCost,
    NodeOrd: NodeOrdStrategy<Cost, PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<Cost>,
    Chaining: ChainingStrategy<Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy<Cost>,
    PrimaryMatch: PrimaryMatchStrategy<Cost>,
    PrimaryRange: PrimaryRangeStrategy,
> Eq
    for AlignmentStrategySelection<
        AlphabetType,
        Cost,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
        PrimaryRange,
    >
{
}
