use std::fmt::Debug;

use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::{
    a_star_aligner::{
        template_switch_distance::strategies::{
            chaining::{ChainingStrategy, NoChainingStrategy, PrecomputeOnlyChainingStrategy},
            node_ord::{AntiDiagonalNodeOrdStrategy, CostOnlyNodeOrdStrategy, NodeOrdStrategy},
            primary_match::AllowPrimaryMatchStrategy,
            secondary_deletion::AllowSecondaryDeletionStrategy,
            shortcut::NoShortcutStrategy,
            template_switch_count::NoTemplateSwitchCountStrategy,
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                TemplateSwitchMinLengthStrategy,
            },
            AlignmentStrategySelection,
        },
        template_switch_distance_a_star_align,
    },
    config::TemplateSwitchConfig,
};
use log::info;

use crate::Cli;

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchNodeOrdStrategySelector {
    CostOnly,
    AntiDiagonal,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchMinLengthStrategySelector {
    None,
    Lookahead,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchChainingStrategySelector {
    None,
    PrecomputeOnly,
}

pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    align_a_star_template_switch_distance_select_node_ord_strategy(cli, reference, query);
}

fn align_a_star_template_switch_distance_select_node_ord_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    match cli.ts_node_ord_strategy {
        TemplateSwitchNodeOrdStrategySelector::CostOnly => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                CostOnlyNodeOrdStrategy,
            >(cli, reference, query)
        }
        TemplateSwitchNodeOrdStrategySelector::AntiDiagonal => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                AntiDiagonalNodeOrdStrategy,
            >(cli, reference, query)
        }
    }
}

fn align_a_star_template_switch_distance_select_template_switch_min_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    match cli.ts_min_length_strategy {
        TemplateSwitchMinLengthStrategySelector::None => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                NoTemplateSwitchMinLengthStrategy,
            >(cli, reference, query)
        }
        TemplateSwitchMinLengthStrategySelector::Lookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                LookaheadTemplateSwitchMinLengthStrategy,
            >(cli, reference, query)
        }
    }
}

fn align_a_star_template_switch_select_chaining_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    match cli.ts_chaining_strategy {
        TemplateSwitchChainingStrategySelector::None => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                NoChainingStrategy,
            >(cli, reference, query)
        }
        TemplateSwitchChainingStrategySelector::PrecomputeOnly => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy,
            >(cli, reference, query)
        }
    }
}

fn align_a_star_template_switch_distance_call<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy,
    Chaining: ChainingStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    let mut config_path = cli.configuration_directory.clone();
    info!("Loading alignment config directory {config_path:?}");

    config_path.push("config.tsa");
    let config_file = std::io::BufReader::new(std::fs::File::open(config_path).unwrap());
    let costs = TemplateSwitchConfig::read_plain(config_file).unwrap();

    info!("Calling aligner...");
    let alignment = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            AlphabetType,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy,
            AllowPrimaryMatchStrategy,
        >,
        _,
    >(reference, query, costs);

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{}", alignment);
}
