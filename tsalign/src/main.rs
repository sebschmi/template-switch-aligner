use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use clap::{Parser, ValueEnum};
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
use template_switch_distance_type_selectors::{
    align_a_star_template_switch_distance, TemplateSwitchNodeOrdStrategy,
};

mod template_switch_distance_type_selectors;

#[derive(Parser)]
struct Cli {
    /// The path to the reference fasta file.
    #[clap(long, short, requires = "query")]
    reference: Option<PathBuf>,

    /// The path to the query fasta file.
    #[clap(long, short, requires = "reference")]
    query: Option<PathBuf>,

    /// The path to a fasta file containing both the reference and the query.
    #[clap(long, short, conflicts_with_all = ["reference", "query"])]
    twin_fasta: Option<PathBuf>,

    /// A directory containing the configuration files.
    ///
    /// See the README for its layout.
    #[clap(long, short, default_value = "sample_tsa_config")]
    configuration_directory: PathBuf,

    #[clap(long, default_value = "a-star-template-switch")]
    alignment_method: AlignmentMethod,

    #[clap(long, default_value = "anti-diagonal")]
    ts_node_ord_strategy: TemplateSwitchNodeOrdStrategy,
}

#[derive(Clone, ValueEnum)]
enum AlignmentMethod {
    Matrix,
    AStarGapAffine,
    AStarTemplateSwitch,
}

fn main() {
    let cli = Cli::parse();

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();

    let (reference, query) = if let Some(twin_fasta) = &cli.twin_fasta {
        let sequences = read_fasta_file(twin_fasta, &mut sequence_store).unwrap();

        assert_eq!(
            sequences.len(),
            2,
            "Twin sequence file contains not exactly two records"
        );
        (
            sequence_store.get(&sequences[0].sequence_handle),
            sequence_store.get(&sequences[1].sequence_handle),
        )
    } else if let (Some(reference), Some(query)) = (&cli.reference, &cli.query) {
        let reference_sequences = read_fasta_file(reference, &mut sequence_store).unwrap();
        let query_sequences = read_fasta_file(query, &mut sequence_store).unwrap();

        assert_eq!(
            reference_sequences.len(),
            1,
            "Reference sequence file contains not exactly one record"
        );
        assert_eq!(
            query_sequences.len(),
            1,
            "Query sequence file contains not exactly one record"
        );
        (
            sequence_store.get(&reference_sequences[0].sequence_handle),
            sequence_store.get(&query_sequences[0].sequence_handle),
        )
    } else {
        panic!("No fasta input file given")
    };

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

    println!("{}", alignment);
}
