use std::path::PathBuf;

use clap::Parser;
use compact_genome::{
    implementation::DefaultSequenceStore,
    interface::{alphabet::dna_alphabet::DnaAlphabet, sequence_store::SequenceStore},
    io::fasta::read_fasta_file,
};
use lib_tsalign::{
    alignment_configuration::AlignmentConfiguration, alignment_matrix::AlignmentMatrix,
};
use traitsequence::interface::Sequence;

#[derive(Parser)]
struct Cli {
    /// The path to the reference fasta file.
    #[clap(long, short)]
    reference: PathBuf,

    /// The path to the query fasta file.
    #[clap(long, short)]
    query: PathBuf,

    #[clap(long, short, default_value = "1")]
    match_score: i64,

    #[clap(long, short, default_value = "-1")]
    substitution_score: i64,

    #[clap(long, short, default_value = "-2")]
    insertion_score: i64,

    #[clap(long, short, default_value = "-2")]
    deletion_score: i64,
}

fn main() {
    let cli = Cli::parse();
    let configuration = AlignmentConfiguration {
        match_score: cli.match_score.into(),
        substitution_score: cli.substitution_score.into(),
        insertion_score: cli.insertion_score.into(),
        deletion_score: cli.deletion_score.into(),
    };

    let mut sequence_store = DefaultSequenceStore::<DnaAlphabet>::new();
    let reference_sequences = read_fasta_file(&cli.reference, &mut sequence_store).unwrap();
    let query_sequences = read_fasta_file(&cli.query, &mut sequence_store).unwrap();

    let reference = sequence_store.get(&reference_sequences[0].sequence_handle);
    let query = sequence_store.get(&query_sequences[0].sequence_handle);

    let mut alignment_matrix = AlignmentMatrix::new(configuration, reference.len(), query.len());
    let score = alignment_matrix.align(reference, query);
    println!("Score: {}", score.as_i64());
}
