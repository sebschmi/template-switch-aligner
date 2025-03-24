use std::{
    fmt::Debug,
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use clap::{Args, Parser, ValueEnum};
use compact_genome::{
    implementation::{
        alphabets::{
            dna_alphabet::DnaAlphabet, dna_alphabet_or_n::DnaAlphabetOrN,
            dna_iupac_nucleic_acid_alphabet::DnaIupacNucleicAcidAlphabet,
            rna_alphabet::RnaAlphabet, rna_alphabet_or_n::RnaAlphabetOrN,
            rna_iupac_nucleic_acid_alphabet::RnaIupacNucleicAcidAlphabet,
        },
        DefaultSequenceStore,
    },
    interface::{alphabet::Alphabet, sequence::GenomeSequence, sequence_store::SequenceStore},
    io::fasta::read_fasta_file,
};
use lib_tsalign::{
    a_star_aligner::{gap_affine_edit_distance, gap_affine_edit_distance_a_star_align},
    alignment_configuration::AlignmentConfiguration,
    alignment_matrix::AlignmentMatrix,
    costs::U64Cost,
};
use log::{debug, info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use template_switch_distance_type_selectors::{
    align_a_star_template_switch_distance, TemplateSwitchChainingStrategySelector,
    TemplateSwitchMinLengthStrategySelector, TemplateSwitchNodeOrdStrategySelector,
};

mod template_switch_distance_type_selectors;

#[derive(Parser)]
pub struct Cli {
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

    /// If set, template switches are not allowed.
    ///
    /// Use this to compare a template switch alignment against an alignment with out template switches.
    #[clap(long)]
    no_ts: bool,

    /// A cost limit for the alignment.
    ///
    /// If there is no alignment with that cost, the aligner will abort without result.
    #[clap(long)]
    cost_limit: Option<U64Cost>,

    /// An approximate memory limit in bytes for the aligner.
    ///
    /// If it is exceeded, then the aligner will abort without result.
    #[clap(long)]
    memory_limit: Option<usize>,
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
    Rna,
    RnaN,
    DnaIupac,
    RnaIupac,
}

pub fn cli(cli: Cli) {
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

    match cli.alphabet {
        InputAlphabet::Dna => execute_with_alphabet::<DnaAlphabet>(cli),
        InputAlphabet::DnaN => execute_with_alphabet::<DnaAlphabetOrN>(cli),
        InputAlphabet::Rna => execute_with_alphabet::<RnaAlphabet>(cli),
        InputAlphabet::RnaN => execute_with_alphabet::<RnaAlphabetOrN>(cli),
        InputAlphabet::DnaIupac => execute_with_alphabet::<DnaIupacNucleicAcidAlphabet>(cli),
        InputAlphabet::RnaIupac => execute_with_alphabet::<RnaIupacNucleicAcidAlphabet>(cli),
    }
}

fn execute_with_alphabet<AlphabetType: Alphabet + Debug + Clone + Eq + 'static>(cli: Cli) {
    let mut skip_characters = Vec::new();
    for character in cli.skip_characters.bytes().map(usize::from) {
        if skip_characters.len() <= character {
            skip_characters.resize(character + 1, false);
        }
        skip_characters[character] = true;
    }
    let skip_characters = skip_characters;

    let mut sequence_store = DefaultSequenceStore::<AlphabetType>::new();
    let sequences = if let Some(pair_fasta) = &cli.input.pair_fasta {
        info!("Loading pair file {pair_fasta:?}");
        let sequences = read_fasta_file(
            pair_fasta,
            &mut sequence_store,
            false,
            true,
            &skip_characters,
        )
        .unwrap_or_else(|error| panic!("Error loading pair file: {error}"));

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
        AlignmentMethod::AStarTemplateSwitch => align_a_star_template_switch_distance(
            cli,
            reference,
            query,
            &format!("{} {}", sequences[0].id, sequences[0].comment),
            &format!("{} {}", sequences[1].id, sequences[1].comment),
        ),
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

    let mut alignment_matrix =
        AlignmentMatrix::<U64Cost>::new(configuration, reference.len(), query.len());
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
        gap_affine_edit_distance::ScoringTable::<U64Cost> {
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
