use clap::Parser;
use compact_genome::{
    implementation::DefaultSequenceStore,
    interface::{alphabet::dna_alphabet::DnaAlphabet, sequence_store::SequenceStore},
};
use lib_tsalign::{
    alignment_configuration::AlignmentConfiguration, alignment_matrix::AlignmentMatrix,
};
use traitsequence::interface::Sequence;

#[derive(Parser)]
struct Cli {
    #[clap(long, short)]
    reference: String,

    #[clap(long, short)]
    query: String,

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
    let reference = sequence_store
        .add_from_slice_u8(cli.reference.as_bytes())
        .unwrap();
    let query = sequence_store
        .add_from_slice_u8(cli.query.as_bytes())
        .unwrap();
    let reference = sequence_store.get(&reference);
    let query = sequence_store.get(&query);

    let mut alignment_matrix = AlignmentMatrix::new(configuration, reference.len(), query.len());
    let score = alignment_matrix.align(reference, query);
    println!("Score: {}", score.as_i64());
}
