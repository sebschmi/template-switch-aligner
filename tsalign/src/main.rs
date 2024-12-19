#![deny(clippy::mod_module_files)]

use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use clap::{Args, Parser, ValueEnum};
use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, DefaultSequenceStore},
    interface::{alphabet::Alphabet, sequence::GenomeSequence, sequence_store::SequenceStore},
    io::fasta::read_fasta_file,
};
use lib_tsalign::{
    a_star_aligner::{gap_affine_edit_distance, gap_affine_edit_distance_a_star_align},
    alignment_configuration::AlignmentConfiguration,
    alignment_matrix::AlignmentMatrix,
};
use log::{debug, info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use template_switch_distance_type_selectors::{
    align_a_star_template_switch_distance, TemplateSwitchChainingStrategySelector,
    TemplateSwitchMinLengthStrategySelector, TemplateSwitchNodeOrdStrategySelector,
};

mod template_switch_distance_type_selectors;

#[derive(Parser)]
struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    #[command(flatten)]
    input: CliInput,

    /// The file to store the alignment statistics in toml format.
    #[clap(long, short = 'o')]
    output: Option<PathBuf>,

    /// The alphabet present in the input files.
    ///
    /// This must also match the alphabet used in the config.
    #[clap(long, short = 'a', default_value = "dna-n")]
    alphabet: InputAlphabet,

    /// A string of (ASCII) characters that should be skipped in the input fasta.
    ///
    /// For example, `-` characters caused by alignment hints can be skipped this way.
    #[clap(long, default_value = "")]
    skip_characters: String,

    /// A directory containing the configuration files.
    ///
    /// See the README for its layout.
    #[clap(long, short = 'c', default_value = "sample_tsa_config")]
    configuration_directory: PathBuf,

    #[clap(long, default_value = "a-star-template-switch")]
    alignment_method: AlignmentMethod,

    #[clap(long, default_value = "anti-diagonal")]
    ts_node_ord_strategy: TemplateSwitchNodeOrdStrategySelector,

    #[clap(long, default_value = "lookahead")]
    ts_min_length_strategy: TemplateSwitchMinLengthStrategySelector,

    #[clap(long, default_value = "none")]
    ts_chaining_strategy: TemplateSwitchChainingStrategySelector,
}

#[derive(Args)]
#[group(required = true)]
struct CliInput {
    /// The path to the reference fasta file.
    #[clap(long, short = 'r', requires = "query", group = "input")]
    reference: Option<PathBuf>,

    /// The path to the query fasta file.
    #[clap(long, short = 'q', requires = "reference", group = "input")]
    query: Option<PathBuf>,

    /// The path to a fasta file containing both the reference and the query.
    #[clap(long, short = 'p', conflicts_with_all = ["reference", "query"], group = "input")]
    pair_fasta: Option<PathBuf>,
}

#[derive(Clone, PartialEq, Eq, ValueEnum)]
enum AlignmentMethod {
    Matrix,
    AStarGapAffine,
    AStarTemplateSwitch,
}

#[derive(Debug, Clone, Eq, PartialEq, ValueEnum)]
enum InputAlphabet {
    Dna,
    DnaN,
}

fn main() {
    let cli = Cli::parse();

    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    if cli.alignment_method != AlignmentMethod::AStarTemplateSwitch
        && cli.alphabet != InputAlphabet::Dna
    {
        panic!("Unsupported alphabet type: {:?}", cli.alphabet);
    }

    let mut skip_characters = Vec::new();
    for character in cli.skip_characters.bytes().map(usize::from) {
        if skip_characters.len() <= character {
            skip_characters.resize(character + 1, false);
        }
        skip_characters[character] = true;
    }
    let skip_characters = skip_characters;

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();
    let sequences = if let Some(pair_fasta) = &cli.input.pair_fasta {
        info!("Loading pair file {pair_fasta:?}");
        let sequences = read_fasta_file(
            pair_fasta,
            &mut sequence_store,
            false,
            true,
            &skip_characters,
        )
        .unwrap();

        assert_eq!(
            sequences.len(),
            2,
            "Pair sequence file contains not exactly two records"
        );

        sequences
    } else if let (Some(reference), Some(query)) = (&cli.input.reference, &cli.input.query) {
        info!("Loading reference file {reference:?}");
        let mut sequences = read_fasta_file(
            reference,
            &mut sequence_store,
            false,
            true,
            &skip_characters,
        )
        .unwrap();
        assert_eq!(
            sequences.len(),
            1,
            "Reference sequence file contains not exactly one record"
        );

        info!("Loading query file {query:?}");
        sequences.extend(
            read_fasta_file(query, &mut sequence_store, false, true, &skip_characters).unwrap(),
        );
        assert_eq!(
            sequences.len(),
            1,
            "Query sequence file contains not exactly one record"
        );

        sequences
    } else {
        panic!("No fasta input file given")
    };

    let reference = sequence_store.get(&sequences[0].sequence_handle);
    let query = sequence_store.get(&sequences[1].sequence_handle);

    debug!("Choosing alignment method...");
    match cli.alignment_method {
        AlignmentMethod::Matrix => align_matrix(cli, reference, query),
        AlignmentMethod::AStarGapAffine => {
            align_a_star_gap_affine_edit_distance(cli, reference, query)
        }
        AlignmentMethod::AStarTemplateSwitch => {
            align_a_star_template_switch_distance(cli, reference, query)
        }
    }
}

fn align_matrix<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    if cli.output.is_some() {
        panic!("Outputting statistics not supported by matrix alignment");
    }

    #[derive(serde::Deserialize)]
    struct MatrixConfig {
        match_cost: u64,
        substitution_cost: u64,
        indel_cost: u64,
    }

    let mut config_path = cli.configuration_directory.clone();
    config_path.push("matrix.toml");
    let mut config_file = BufReader::new(File::open(config_path).unwrap());
    let mut config = String::new();
    config_file.read_to_string(&mut config).unwrap();
    let matrix_config: MatrixConfig = toml::from_str(&config).unwrap();

    let configuration = AlignmentConfiguration {
        match_cost: matrix_config.match_cost.into(),
        substitution_cost: matrix_config.substitution_cost.into(),
        insertion_cost: matrix_config.indel_cost.into(),
        deletion_cost: matrix_config.indel_cost.into(),
    };

    let mut alignment_matrix = AlignmentMatrix::new(configuration, reference.len(), query.len());
    let cost = alignment_matrix.align(reference, query);
    println!("Cost: {}", cost);
}

fn align_a_star_gap_affine_edit_distance<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    #[derive(serde::Deserialize)]
    struct GapAffineConfig {
        match_cost: u64,
        substitution_cost: u64,
        gap_open_cost: u64,
        gap_extend_cost: u64,
    }

    let mut config_path = cli.configuration_directory.clone();
    config_path.push("a_star_gap_affine.toml");
    let mut config_file = BufReader::new(File::open(config_path).unwrap());
    let mut config = String::new();
    config_file.read_to_string(&mut config).unwrap();
    let gap_affine_config: GapAffineConfig = toml::from_str(&config).unwrap();

    let alignment = gap_affine_edit_distance_a_star_align(
        reference,
        query,
        gap_affine_edit_distance::ScoringTable {
            match_cost: gap_affine_config.match_cost.into(),
            substitution_cost: gap_affine_config.substitution_cost.into(),
            gap_open_cost: gap_affine_config.gap_open_cost.into(),
            gap_extend_cost: gap_affine_config.gap_extend_cost.into(),
        },
    );

    if let Some(output) = cli.output {
        use std::io::Write;
        let mut output = std::io::BufWriter::new(std::fs::File::create(output).unwrap());
        write!(output, "{}", toml::to_string(&alignment).unwrap()).unwrap();
    }

    println!("{}", alignment);
}
