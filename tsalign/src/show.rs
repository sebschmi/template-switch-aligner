use std::{fs::File, io::Read, path::PathBuf};

use clap::Parser;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{a_star_sequences::SequencePair, AlignmentResult},
        template_switch_distance::AlignmentType,
    },
    costs::U64Cost,
};
use log::{info, warn, LevelFilter};
use parse_template_switches::TSShow;
use simplelog::{ColorChoice, TermLogger, TerminalMode};

mod alignment_stream;
mod parse_template_switches;

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

fn show_template_switches(result: &AlignmentResult<AlignmentType, U64Cost>) {
    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        warn!("Alignment was aborted early, no template switches present");
        return;
    };

    info!("Collecting template switches");
    let template_switches = parse_template_switches::parse(alignment);
    info!("Found {} template switches", template_switches.len());

    for (index, template_switch) in template_switches.iter().enumerate() {
        info!("Showing template switch {}", index + 1);
        show_template_switch(template_switch, &statistics.sequences);
    }
}

fn show_template_switch(_template_switch: &TSShow<AlignmentType>, _sequences: &SequencePair) {
    todo!()
}
