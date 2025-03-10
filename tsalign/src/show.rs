use std::{fs::File, io::Read, path::PathBuf};

use clap::Parser;
use lib_tsalign::{
    a_star_aligner::{alignment_result::AlignmentResult, gap_affine_edit_distance::AlignmentType},
    costs::U64Cost,
};
use log::{info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// Path to a toml output file of tsalign.
    #[clap(long, short = 'i')]
    input: PathBuf,
}

pub fn cli(cli: Cli) {
    TermLogger::init(
        cli.log_level,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

    info!("Reading tsalign output toml file {:?}", cli.input);
    let mut input = String::new();
    File::open(cli.input)
        .unwrap_or_else(|error| panic!("Error opening input file: {error}"))
        .read_to_string(&mut input)
        .unwrap_or_else(|error| panic!("Error reading input file: {error}"));

    let result =
        toml::from_str(&input).unwrap_or_else(|error| panic!("Error parsing input file: {error}"));
    show_template_switches(&result);
}

fn show_template_switches(_result: &AlignmentResult<AlignmentType, U64Cost>) {
    todo!()
}
