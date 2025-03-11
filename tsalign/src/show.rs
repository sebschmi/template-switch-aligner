use std::{
    fs::File,
    io::{stdout, Read},
    path::PathBuf,
};

use alignment_stream::AlignmentCoordinates;
use clap::Parser;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{a_star_sequences::SequencePair, AlignmentResult},
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{info, warn, LevelFilter};
use mutlipair_alignment_renderer::MultipairAlignmentRenderer;
use parse_template_switches::TSShow;
use simplelog::{ColorChoice, TermLogger, TerminalMode};

mod alignment_stream;
mod mutlipair_alignment_renderer;
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

fn show_template_switch(template_switch: &TSShow<AlignmentType>, sequences: &SequencePair) {
    // println!("Showing template switch\n{template_switch:?}");

    let reference = &sequences.reference;
    let reference_c: String = sequences.reference_rc.chars().rev().collect();
    let query = &sequences.query;
    let query_c: String = sequences.query_rc.chars().rev().collect();

    let (
        primary_name,
        primary,
        primary_c,
        primary_coordinate_picker,
        anti_primary_name,
        anti_primary,
        anti_primary_c,
        anti_primary_coordinate_picker,
        invert_alignment,
    ) = match template_switch.primary {
        TemplateSwitchPrimary::Reference => (
            "Ref".to_string(),
            reference,
            &reference_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| {
                alignment_coordinates.reference()
            }) as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            "Qry".to_string(),
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| alignment_coordinates.query())
                as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            false,
        ),
        TemplateSwitchPrimary::Query => (
            "Qry".to_string(),
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| alignment_coordinates.query())
                as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            "Ref".to_string(),
            reference,
            &reference_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| {
                alignment_coordinates.reference()
            }) as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            true,
        ),
    };
    let primary_equals_secondary = (template_switch.primary == TemplateSwitchPrimary::Reference
        && template_switch.secondary == TemplateSwitchSecondary::Reference)
        || (template_switch.primary == TemplateSwitchPrimary::Query
            && template_switch.secondary == TemplateSwitchSecondary::Query);

    if primary_equals_secondary {
        todo!()
    } else {
        let primary_forward_name = format!("{primary_name}F");
        let primary_reverse_name = format!("{primary_name}R");
        let f1_name = format!("{anti_primary_name}1");
        let f2_name = format!("{anti_primary_name}2");
        let f3_name = format!("{anti_primary_name}3");

        let primary_offset = primary_coordinate_picker(&template_switch.upstream_offset);
        let primary_limit = primary_coordinate_picker(&template_switch.downstream_limit);
        let anti_primary_offset = anti_primary_coordinate_picker(&template_switch.upstream_offset);
        let anti_primary_limit = anti_primary_coordinate_picker(&template_switch.downstream_limit);

        let ts_reverse: String =
            anti_primary[anti_primary_coordinate_picker(&template_switch.sp1_offset)
                ..anti_primary_coordinate_picker(&template_switch.sp4_offset)]
                .chars()
                .rev()
                .collect();
        let ts_reverse_alignment: Vec<_> = template_switch
            .template_switch
            .iter()
            .copied()
            .rev()
            .collect();

        let mut renderer = MultipairAlignmentRenderer::new(
            primary_forward_name.clone(),
            &primary[primary_offset..primary_limit],
        );
        renderer.add_aligned_sequence(
            &primary_forward_name,
            0,
            primary_reverse_name.clone(),
            &primary_c[primary_offset..primary_limit],
            &[(primary_limit - primary_offset, AlignmentType::PrimaryMatch)],
            false,
            false,
        );

        renderer.add_aligned_sequence(
            &primary_forward_name,
            0,
            f1_name.clone(),
            &anti_primary
                [anti_primary_offset..anti_primary_coordinate_picker(&template_switch.sp1_offset)],
            &template_switch.upstream,
            true,
            invert_alignment,
        );
        renderer.add_aligned_sequence(
            &primary_forward_name,
            primary_coordinate_picker(&template_switch.sp4_offset) - primary_offset,
            f3_name.clone(),
            &anti_primary
                [anti_primary_coordinate_picker(&template_switch.sp4_offset)..anti_primary_limit],
            &template_switch.downstream,
            true,
            invert_alignment,
        );

        renderer.add_aligned_sequence(
            &primary_reverse_name,
            template_switch.sp3_primary_offset - primary_offset,
            f2_name.clone(),
            &ts_reverse,
            &ts_reverse_alignment,
            true,
            invert_alignment,
        );

        renderer.render(stdout()).unwrap();
    }
}
