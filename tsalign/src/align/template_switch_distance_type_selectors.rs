use std::fmt::Debug;

use clap::ValueEnum;
use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use lib_tsalign::{
    a_star_aligner::{
        alignment_geometry::{AlignmentCoordinates, AlignmentRange},
        template_switch_distance::strategies::{
            AlignmentStrategySelection,
            chaining::{
                ChainingStrategy, LowerBoundChainingStrategy, NoChainingStrategy,
                PrecomputeOnlyChainingStrategy,
            },
            node_ord::{AntiDiagonalNodeOrdStrategy, CostOnlyNodeOrdStrategy, NodeOrdStrategy},
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
            template_switch_total_length::MaxTemplateSwitchTotalLengthStrategy,
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
    LowerBound,
}

pub fn align_a_star_template_switch_distance<
    AlphabetType: Alphabet + Debug + Clone + Eq,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
    reference_name: &str,
    query_name: &str,
) {
    align_a_star_template_switch_distance_select_node_ord_strategy(
        cli,
        reference,
        query,
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
    reference_name: &str,
    query_name: &str,
) {
    match cli.ts_node_ord_strategy {
        TemplateSwitchNodeOrdStrategySelector::CostOnly => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                CostOnlyNodeOrdStrategy,
            >(cli, reference, query, reference_name, query_name)
        }
        TemplateSwitchNodeOrdStrategySelector::AntiDiagonal => {
            align_a_star_template_switch_distance_select_template_switch_min_length_strategy::<
                _,
                _,
                AntiDiagonalNodeOrdStrategy,
            >(cli, reference, query, reference_name, query_name)
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
            >(cli, reference, query, reference_name, query_name)
        }
        TemplateSwitchMinLengthStrategySelector::Lookahead => {
            align_a_star_template_switch_select_chaining_strategy::<
                _,
                _,
                NodeOrd,
                LookaheadTemplateSwitchMinLengthStrategy<U64Cost>,
            >(cli, reference, query, reference_name, query_name)
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
            >(cli, reference, query, reference_name, query_name)
        }
        TemplateSwitchChainingStrategySelector::PrecomputeOnly => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                PrecomputeOnlyChainingStrategy<U64Cost>,
            >(cli, reference, query, reference_name, query_name)
        }
        TemplateSwitchChainingStrategySelector::LowerBound => {
            align_a_star_template_switch_select_no_ts_strategy::<
                _,
                _,
                NodeOrd,
                TemplateSwitchMinLength,
                LowerBoundChainingStrategy<U64Cost>,
            >(cli, reference, query, reference_name, query_name)
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
    reference_name: &str,
    query_name: &str,
) {
    if cli.no_ts {
        align_a_star_template_switch_distance_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            MaxTemplateSwitchCountStrategy,
        >(cli, reference, query, reference_name, query_name, 0)
    } else {
        align_a_star_template_switch_distance_call::<
            _,
            _,
            NodeOrd,
            TemplateSwitchMinLength,
            Chaining,
            NoTemplateSwitchCountStrategy,
        >(cli, reference, query, reference_name, query_name, ())
    }
}

fn align_a_star_template_switch_distance_call<
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

    let range = Some(parse_range(&cli, reference.len(), query.len()));

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
            MaxTemplateSwitchTotalLengthStrategy,
        >,
        _,
    >(
        reference,
        query,
        reference_name,
        query_name,
        range,
        costs,
        cli.cost_limit,
        cli.memory_limit,
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

fn parse_range(cli: &Cli, reference_length: usize, query_length: usize) -> AlignmentRange {
    let complete_reference_range = 0..reference_length;
    let complete_query_range = 0..query_length;

    let (reference_range, query_range) = if let Some(rq_ranges) = cli.rq_ranges.as_ref() {
        let mut rq_ranges = rq_ranges.chars().peekable();

        let mut reference_range = None;
        let mut query_range = None;

        while rq_ranges.peek().is_some() {
            let rq = rq_ranges.next().unwrap();

            while let Some(c) = rq_ranges.peek() {
                if c.is_whitespace() {
                    rq_ranges.next().unwrap();
                } else {
                    break;
                }
            }

            let mut offset = String::new();
            while let Some(c) = rq_ranges.peek() {
                if c.is_numeric() {
                    offset.push(rq_ranges.next().unwrap());
                } else {
                    break;
                }
            }

            // Parse ..
            assert_eq!(rq_ranges.next(), Some('.'));
            assert_eq!(rq_ranges.next(), Some('.'));

            let mut limit = String::new();
            while let Some(c) = rq_ranges.peek() {
                if c.is_numeric() {
                    limit.push(rq_ranges.next().unwrap());
                } else {
                    break;
                }
            }

            let offset = offset.parse().unwrap();
            let limit = limit.parse().unwrap();

            match rq {
                'R' => {
                    assert!(reference_range.is_none());
                    reference_range = Some(offset..limit)
                }
                'Q' => {
                    assert!(query_range.is_none());
                    query_range = Some(offset..limit)
                }
                _ => panic!(),
            }
        }

        assert!(
            reference_range.is_none()
                || (cli.reference_offset.is_none() && cli.reference_limit.is_none())
        );
        assert!(query_range.is_none() || (cli.query_offset.is_none() && cli.query_limit.is_none()));

        (
            reference_range.unwrap_or(complete_reference_range),
            query_range.unwrap_or(complete_query_range),
        )
    } else {
        (complete_reference_range, complete_query_range)
    };

    AlignmentRange::new_offset_limit(
        AlignmentCoordinates::new(
            cli.reference_offset.unwrap_or(reference_range.start),
            cli.query_offset.unwrap_or(query_range.start),
        ),
        AlignmentCoordinates::new(
            cli.reference_limit.unwrap_or(reference_range.end),
            cli.query_limit.unwrap_or(query_range.end),
        ),
    )
}
