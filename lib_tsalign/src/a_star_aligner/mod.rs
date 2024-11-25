use std::{fmt::Debug, time::Instant};

use alignment_result::{AlignmentResult, IAlignmentType};
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use generic_a_star::{AStar, AStarContext, AStarNode, AStarResult};
use template_switch_distance::strategies::{
    template_switch_count::NoTemplateSwitchCountStrategy, AlignmentStrategySelector,
};
use traitsequence::interface::Sequence;

use crate::config;

pub mod alignment_result;
pub mod gap_affine_edit_distance;
pub mod template_switch_distance;
#[cfg(test)]
mod tests;

pub trait AlignmentContext: AStarContext {
    type AlphabetType: Alphabet;

    type AlignmentType: IAlignmentType
        + Debug
        + From<<<Self as AStarContext>::Node as AStarNode>::EdgeType>;

    type SubsequenceType: GenomeSequence<Self::AlphabetType, Self::SubsequenceType> + ?Sized;

    fn reference(&self) -> &Self::SubsequenceType;

    fn query(&self) -> &Self::SubsequenceType;
}

fn a_star_align<Context: AStarContext + AlignmentContext>(
    context: Context,
) -> AlignmentResult<Context::AlignmentType>
where
    <Context::Node as AStarNode>::EdgeType: IAlignmentType,
{
    let start_time = Instant::now();

    // Perform forwards search.
    let mut a_star = AStar::new(context);
    a_star.initialise();
    let AStarResult::FoundTarget { cost, .. } = a_star.search() else {
        unreachable!("The search can always reach a target.")
    };

    let mut alignment = Vec::new();

    // Backtrack.
    for alignment_type in a_star
        .backtrack()
        .into_iter()
        .map(<Context as AlignmentContext>::AlignmentType::from)
    {
        if !alignment_type.is_internal() {
            if let Some((count, previous_alignment_type)) = alignment.last_mut() {
                if alignment_type.is_repeated(previous_alignment_type) {
                    *count += 1;
                } else {
                    alignment.push((1, alignment_type));
                }
            } else {
                alignment.push((1, alignment_type));
            }
        }
    }

    alignment.reverse();

    let end_time = Instant::now();
    let duration = (end_time - start_time).as_secs_f64();

    AlignmentResult::new(
        alignment,
        cost,
        duration,
        a_star.performance_counters().opened_nodes,
        a_star.performance_counters().closed_nodes,
        a_star.performance_counters().suboptimal_opened_nodes,
        a_star.context().reference().len(),
        a_star.context().query().len(),
    )
}

pub fn gap_affine_edit_distance_a_star_align<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    scoring_table: gap_affine_edit_distance::ScoringTable,
) -> AlignmentResult<gap_affine_edit_distance::AlignmentType> {
    a_star_align(gap_affine_edit_distance::Context::new(
        reference,
        query,
        scoring_table,
    ))
}

pub fn template_switch_distance_a_star_align<
    Strategies: AlignmentStrategySelector<TemplateSwitchCount = NoTemplateSwitchCountStrategy>,
    SubsequenceType: GenomeSequence<Strategies::Alphabet, SubsequenceType> + ?Sized,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    config: config::TemplateSwitchConfig<Strategies::Alphabet>,
) -> AlignmentResult<template_switch_distance::AlignmentType> {
    a_star_align(template_switch_distance::Context::<
        SubsequenceType,
        Strategies,
    >::new(reference, query, config, ()))
}
