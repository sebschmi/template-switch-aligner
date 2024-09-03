use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use compact_genome::{
    implementation::DefaultSequenceStore,
    interface::{
        alphabet::{dna_alphabet::DnaAlphabet, Alphabet},
        sequence::GenomeSequence,
        sequence_store::SequenceStore,
    },
    io::fasta::read_fasta_file,
};
use lib_tsalign::{
    a_star_aligner::{
        gap_affine_edit_distance::ScoringTable, gap_affine_edit_distance_a_star_align,
    },
    alignment_configuration::AlignmentConfiguration,
    alignment_matrix::AlignmentMatrix,
};

#[derive(Parser)]
struct Cli {
    /// The path to the reference fasta file.
    #[clap(long, short)]
    reference: PathBuf,

    /// The path to the query fasta file.
    #[clap(long, short)]
    query: PathBuf,

    #[clap(long, short, allow_negative_numbers = true, default_value = "1")]
    match_score: i64,

    #[clap(long, short, allow_negative_numbers = true, default_value = "-1")]
    substitution_score: i64,

    #[clap(long, short, allow_negative_numbers = true, default_value = "-2")]
    insertion_score: i64,

    #[clap(long, short, allow_negative_numbers = true, default_value = "-2")]
    deletion_score: i64,

    #[clap(long, short = 'o', allow_negative_numbers = true, default_value = "-3")]
    gap_open_score: i64,

    #[clap(long, short = 'e', allow_negative_numbers = true, default_value = "-1")]
    gap_extend_score: i64,

    #[clap(long, default_value = "a-star-gap-affine-edit-distance")]
    alignment_method: AlignmentMethod,
}

#[derive(Clone, ValueEnum)]
enum AlignmentMethod {
    Matrix,
    AStarGapAffineEditDistance,
}

fn main() {
    let cli = Cli::parse();

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();
    let reference_sequences = read_fasta_file(&cli.reference, &mut sequence_store).unwrap();
    let query_sequences = read_fasta_file(&cli.query, &mut sequence_store).unwrap();

    let reference = sequence_store.get(&reference_sequences[0].sequence_handle);
    let query = sequence_store.get(&query_sequences[0].sequence_handle);

    match cli.alignment_method {
        AlignmentMethod::Matrix => align_matrix(cli, reference, query),
        AlignmentMethod::AStarGapAffineEditDistance => {
            align_a_star_gap_affine_edit_distance(cli, reference, query)
        }
    }
}

fn align_a_star_gap_affine_edit_distance<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    let alignment = gap_affine_edit_distance_a_star_align(
        reference,
        query,
        ScoringTable {
            match_score: cli.match_score.into(),
            substitution_score: cli.substitution_score.into(),
            gap_open_score: cli.gap_open_score.into(),
            gap_extend_score: cli.gap_extend_score.into(),
        },
    );

    println!("Score: {}", alignment.score.as_i64());
    println!(
        "Score per base: {:.2}",
        (alignment.score.as_i64() * 2) as f64 / (reference.len() + query.len()) as f64
    );
    println!("Opened nodes: {}", alignment.opened_nodes);
    println!("Closed nodes: {}", alignment.closed_nodes);
    println!("CIGAR: {}", alignment.cigar());
}

fn align_matrix<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    cli: Cli,
    reference: &SubsequenceType,
    query: &SubsequenceType,
) {
    let configuration = AlignmentConfiguration {
        match_score: cli.match_score.into(),
        substitution_score: cli.substitution_score.into(),
        insertion_score: cli.insertion_score.into(),
        deletion_score: cli.deletion_score.into(),
    };

    let mut alignment_matrix = AlignmentMatrix::new(configuration, reference.len(), query.len());
    let score = alignment_matrix.align(reference, query);
    println!("Score: {}", score.as_i64());
}
