use std::{fmt::Debug, hash::Hash, marker::PhantomData};

use chaining::ChainingStrategy;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use node_ord::NodeOrdStrategy;
use primary_match::PrimaryMatchStrategy;
use secondary_deletion::SecondaryDeletionStrategy;
use shortcut::ShortcutStrategy;
use template_switch_count::TemplateSwitchCountStrategy;
use template_switch_min_length::TemplateSwitchMinLengthStrategy;

use super::{AlignmentType, Context, Identifier};

pub mod chaining;
pub mod node_ord;
pub mod primary_match;
pub mod secondary_deletion;
pub mod shortcut;
pub mod template_switch_count;
pub mod template_switch_min_length;

pub trait AlignmentStrategySelector: Eq + Copy + Ord + std::fmt::Debug {
    type Alphabet: Alphabet;
    type NodeOrd: NodeOrdStrategy<Self::PrimaryMatch>;
    type TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy;
    type Chaining: ChainingStrategy;
    type TemplateSwitchCount: TemplateSwitchCountStrategy;
    type SecondaryDeletion: SecondaryDeletionStrategy;
    type Shortcut: ShortcutStrategy;
    type PrimaryMatch: PrimaryMatchStrategy;
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AlignmentStrategiesNodeMemory<Selector: AlignmentStrategySelector> {
    pub node_ord_strategy: Selector::NodeOrd,
    pub template_switch_min_length_strategy: Selector::TemplateSwitchMinLength,
    pub template_switch_count: Selector::TemplateSwitchCount,
    pub primary_match: Selector::PrimaryMatch,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentStrategiesNodeIdentifier<Strategies: AlignmentStrategySelector> {
    pub primary_match: <Strategies::PrimaryMatch as PrimaryMatchStrategy>::Identifier,
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
        identifier: Identifier<Strategies>,
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
        identifier: Identifier<Strategies>,
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

impl<Strategies: AlignmentStrategySelector> AlignmentStrategiesNodeIdentifier<Strategies> {
    pub fn create_root<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            primary_match: Strategies::PrimaryMatch::create_root_identifier(context),
        }
    }

    pub fn generate_successor<
        SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
    >(
        &self,
        identifier: Identifier<Strategies>,
        alignment_type: AlignmentType,
        context: &Context<'_, '_, SubsequenceType, Strategies>,
    ) -> Self {
        Self {
            primary_match: Strategies::PrimaryMatch::generate_successor_identifier(
                identifier.strategies.primary_match,
                alignment_type,
                context,
            ),
        }
    }
}

pub struct AlignmentStrategySelection<
    AlphabetType: Alphabet,
    NodeOrd: NodeOrdStrategy<PrimaryMatch>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
    Chaining: ChainingStrategy,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    SecondaryDeletion: SecondaryDeletionStrategy,
    Shortcut: ShortcutStrategy,
    PrimaryMatch: PrimaryMatchStrategy,
> {
    #[allow(clippy::type_complexity)]
    phantom_data: PhantomData<(
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    )>,
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > AlignmentStrategySelector
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
    type Alphabet = AlphabetType;
    type NodeOrd = NodeOrd;
    type TemplateSwitchMinLength = TemplateSwitchMinLength;
    type Chaining = Chaining;
    type TemplateSwitchCount = TemplateSwitchCount;
    type SecondaryDeletion = SecondaryDeletion;
    type Shortcut = Shortcut;
    type PrimaryMatch = PrimaryMatch;
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > Debug
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentStrategySelection").finish()
    }
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > Clone
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
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
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > Copy
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > PartialEq
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
    fn eq(&self, other: &Self) -> bool {
        self.phantom_data == other.phantom_data
    }
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > Eq
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > PartialOrd
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.phantom_data.partial_cmp(&other.phantom_data)
    }
}

impl<
        AlphabetType: Alphabet,
        NodeOrd: NodeOrdStrategy<PrimaryMatch>,
        TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
        Chaining: ChainingStrategy,
        TemplateSwitchCount: TemplateSwitchCountStrategy,
        SecondaryDeletion: SecondaryDeletionStrategy,
        Shortcut: ShortcutStrategy,
        PrimaryMatch: PrimaryMatchStrategy,
    > Ord
    for AlignmentStrategySelection<
        AlphabetType,
        NodeOrd,
        TemplateSwitchMinLength,
        Chaining,
        TemplateSwitchCount,
        SecondaryDeletion,
        Shortcut,
        PrimaryMatch,
    >
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.phantom_data.cmp(&other.phantom_data)
    }
}

impl<Strategies: AlignmentStrategySelector> PartialOrd
    for AlignmentStrategiesNodeIdentifier<Strategies>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<Strategies: AlignmentStrategySelector> Ord for AlignmentStrategiesNodeIdentifier<Strategies> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.primary_match.cmp(&other.primary_match)
    }
}

impl<Strategies: AlignmentStrategySelector> Hash for AlignmentStrategiesNodeIdentifier<Strategies> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.primary_match.hash(state);
    }
}

impl<Strategies: AlignmentStrategySelector> Copy for AlignmentStrategiesNodeIdentifier<Strategies> {}
