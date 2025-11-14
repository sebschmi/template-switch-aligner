use std::fmt::Debug;

use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::{
    a_star_aligner::{
        alignment_geometry::AlignmentRange,
        template_switch_distance::strategies::{
            AlignmentStrategySelection,
            chaining::{ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy},
            node_ord::{AntiDiagonalNodeOrdStrategy, NodeOrdStrategy},
            primary_match::AllowPrimaryMatchStrategy,
            primary_range::NoPrunePrimaryRangeStrategy,
            secondary_deletion::AllowSecondaryDeletionStrategy,
            shortcut::NoShortcutStrategy,
            template_switch_count::{
                MaxTemplateSwitchCountStrategy, NoTemplateSwitchCountStrategy,
                TemplateSwitchCountStrategy,
            },
            template_switch_min_length::{
                LookaheadTemplateSwitchMinLengthStrategy, NoTemplateSwitchMinLengthStrategy,
                TemplateSwitchMinLengthStrategy,
            },
            template_switch_total_length::{
                MaxTemplateSwitchTotalLengthStrategy, NoTemplateSwitchTotalLengthStrategy,
                TemplateSwitchTotalLengthStrategy,
            },
        },
        template_switch_distance_a_star_align,
    },
    config::TemplateSwitchConfig,
    costs::U64Cost,
};
use log::info;

use super::Cli;

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchNodeOrdStrategySelector {
    // CostOnly,
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
    // PrecomputeOnly,
    LowerBound,
}

#[derive(Clone, ValueEnum)]
pub enum TemplateSwitchTotalLengthStrategySelector {
    None,
    Maximise,
}

pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
) {
    align_a_star_template_switch_distance_select_node_ord_strategy(
        cli,
        reference,
        query,
        range,
        reference_name,
        query_name,
    );
}

fn align_a_star_template_switch_distance_select_node_ord_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_node_ord_strategy {
        /*TemplateSwitchNodeOrdStrategySelector::CostOnly => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                CostOnlyNodeOrdStrategy,
            >(cli, reference, query, reference_name, query_name)
            unimplemented!("The other option appears to always be better.");
        }*/
        TemplateSwitchNodeOrdStrategySelector::AntiDiagonal => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                AntiDiagonalNodeOrdStrategy,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_distance_select_template_switch_min_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_min_length_strategy {
        TemplateSwitchMinLengthStrategySelector::None => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                NoTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::Lookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_select_chaining_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_chaining_strategy {
        TemplateSwitchChainingStrategySelector::None => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                NoChainingStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
        /*TemplateSwitchChainingStrategySelector::PrecomputeOnly => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy<U64Cost>,
            >(cli, reference, query, reference_name, query_name)
            unimplemented!("No reason to precompute without using the information.");
        }*/
        TemplateSwitchChainingStrategySelector::LowerBound => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                LowerBoundChainingStrategy<U64Cost>,
            >(cli, reference, query, range, reference_name, query_name)
        }
    }
}

fn align_a_star_template_switch_select_no_ts_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
) {
    if cli.no_ts {
        align_a_star_template_switch_select_template_switch_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            MaxTemplateSwitchCountStrategy,
        >(cli, reference, query, range, reference_name, query_name, 0)
    } else {
        align_a_star_template_switch_select_template_switch_total_length_strategy::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
        >(cli, reference, query, range, reference_name, query_name, ())
    }
}

fn align_a_star_template_switch_select_template_switch_total_length_strategy<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) {
    match cli.ts_total_length_strategy {
        TemplateSwitchTotalLengthStrategySelector::None => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                NoTemplateSwitchTotalLengthStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
        TemplateSwitchTotalLengthStrategySelector::Maximise => {
            align_a_star_template_switch_distance_call::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                Chaining,
                TemplateSwitchCount,
                MaxTemplateSwitchTotalLengthStrategy,
            >(
                cli,
                reference,
                query,
                range,
                reference_name,
                query_name,
                template_switch_count_memory,
            )
        }
    }
}

fn align_a_star_template_switch_distance_call<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    NodeOrd: NodeOrdStrategy<U64Cost, AllowPrimaryMatchStrategy>,
    TemplateSwitchMinLength: TemplateSwitchMinLengthStrategy<U64Cost>,
    Chaining: ChainingStrategy<U64Cost>,
    TemplateSwitchCount: TemplateSwitchCountStrategy,
    TemplateSwitchTotalLength: TemplateSwitchTotalLengthStrategy,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    range: Option<AlignmentRange>,
    reference_name: &str,
    query_name: &str,
    template_switch_count_memory: <TemplateSwitchCount as TemplateSwitchCountStrategy>::Memory,
) {
    let mut config_path = cli.configuration_directory.clone();
    info!("Loading alignment config directory {config_path:?}");

    config_path.push("config.tsa");
    let config_file = std::io::BufReader::new(
        std::fs::File::open(&config_path)
            .unwrap_or_else(|error| panic!("Error opening config file {config_path:?}: {error}")),
    );
    let costs = TemplateSwitchConfig::read_plain(config_file)
        .unwrap_or_else(|error| panic!("Error parsing template switch config:\n{error}"));

    info!("Calling aligner...");
    let alignment = template_switch_distance_a_star_align::<
        AlignmentStrategySelection<
            AlphabetType,
            U64Cost,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            TemplateSwitchCount,
            AllowSecondaryDeletionStrategy,
            NoShortcutStrategy<U64Cost>,
            AllowPrimaryMatchStrategy,
            NoPrunePrimaryRangeStrategy,
            TemplateSwitchTotalLength,
        >,
        _,
    >(
        reference,
        query,
        reference_name,
        query_name,
        range,
        &costs,
        cli.cost_limit,
        cli.memory_limit,
        cli.force_label_correcting,
        template_switch_count_memory,
    );
    info!("Finished aligning");

    if let Some(output) = cli.output {
        info!("Outputting alignment statistics to {output:?}");
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{alignment}");
}
