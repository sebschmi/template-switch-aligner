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
    #[clap(long, short, requires = "query")]
    reference: Option<PathBuf>,

    /// The path to the query fasta file.
    #[clap(long, short, requires = "reference")]
    query: Option<PathBuf>,

    /// The path to a fasta file containing both the reference and the query.
    #[clap(long, short, conflicts_with_all = ["reference", "query"])]
    twin_fasta: Option<PathBuf>,

    #[clap(long, short, default_value = "0")]
    match_cost: u64,

    #[clap(long, short, default_value = "2")]
    substitution_cost: u64,

    #[clap(long, short, default_value = "3")]
    insertion_cost: u64,

    #[clap(long, short, default_value = "3")]
    deletion_cost: u64,

    #[clap(long, short = 'o', default_value = "3")]
    gap_open_cost: u64,

    #[clap(long, short = 'e', default_value = "1")]
    gap_extend_cost: u64,

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
            match_cost: cli.match_cost.into(),
            substitution_cost: cli.substitution_cost.into(),
            gap_open_cost: cli.gap_open_cost.into(),
            gap_extend_cost: cli.gap_extend_cost.into(),
        },
    );

    println!("{}", alignment);
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
        match_cost: cli.match_cost.into(),
        substitution_cost: cli.substitution_cost.into(),
        insertion_cost: cli.insertion_cost.into(),
        deletion_cost: cli.deletion_cost.into(),
    };

    let mut alignment_matrix = AlignmentMatrix::new(configuration, reference.len(), query.len());
    let cost = alignment_matrix.align(reference, query);
    println!("Cost: {}", cost);
}
