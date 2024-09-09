use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::a_star_aligner::{
    template_switch_distance::{
        strategies::{
            node_ord::{AntiDiagonalNodeOrdStrategy, CostOnlyNodeOrdStrategy, NodeOrdStrategy},
            AlignmentStrategySelection,
        },
        ScoringTable,
    },
    template_switch_distance_a_star_align,
};

use crate::Cli;

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchNodeOrdStrategy {
    CostOnly,
    AntiDiagonal,
}

pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    align_a_star_template_switch_distance_select_node_ord_strategy(cli, reference, query);
}

fn align_a_star_template_switch_distance_select_node_ord_strategy<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    match cli.ts_node_ord_strategy {
        TemplateSwitchNodeOrdStrategy::CostOnly => {
            align_a_star_template_switch_distance_call::<_, _, CostOnlyNodeOrdStrategy>(
                cli, reference, query,
            )
        }
        TemplateSwitchNodeOrdStrategy::AntiDiagonal => {
            align_a_star_template_switch_distance_call::<_, _, AntiDiagonalNodeOrdStrategy>(
                cli, reference, query,
            )
        }
    }
}

fn align_a_star_template_switch_distance_call<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    let alignment =
        template_switch_distance_a_star_align::<_, _, AlignmentStrategySelection<NodeOrd>>(
            reference,
            query,
            ScoringTable {
                match_cost: cli.match_cost.into(),
                substitution_cost: cli.substitution_cost.into(),
                gap_open_cost: cli.gap_open_cost.into(),
                gap_extend_cost: cli.gap_extend_cost.into(),
            },
        );

    println!("{}", alignment);
}
