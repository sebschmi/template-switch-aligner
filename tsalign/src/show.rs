use std::{
    fs::File,
    io::{Read, stdout},
    path::PathBuf,
};

use alignment_stream::AlignmentCoordinates;
use clap::Parser;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{AlignmentResult, a_star_sequences::SequencePair},
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{LevelFilter, debug, info, warn};
use mutlipair_alignment_renderer::MultipairAlignmentRenderer;
use parse_template_switches::{STREAM_PADDING, TSShow};
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

    debug!("CIGAR: {}", result.cigar());

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
        _primary_c,
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
            true,
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
            false,
        ),
    };
    let primary_equals_secondary = (template_switch.primary == TemplateSwitchPrimary::Reference
        && template_switch.secondary == TemplateSwitchSecondary::Reference)
        || (template_switch.primary == TemplateSwitchPrimary::Query
            && template_switch.secondary == TemplateSwitchSecondary::Query);

    if primary_equals_secondary {
        todo!()
    } else {
        let anti_primary_forward_name = format!("{anti_primary_name}F");
        let anti_primary_reverse_name = format!("{anti_primary_name}R");
        let f1_name = format!("{primary_name}1");
        let f2_name = format!("{primary_name}2");
        let f3_name = format!("{primary_name}3");

        let primary_offset = primary_coordinate_picker(&template_switch.upstream_offset);

        let primary_limit = primary_coordinate_picker(&template_switch.downstream_limit);

        let anti_primary_offset = anti_primary_coordinate_picker(&template_switch.upstream_offset);
        let anti_primary_extended_offset = anti_primary_offset.min(
            template_switch
                .sp3_secondary_offset
                .saturating_sub(STREAM_PADDING),
        );
        let anti_primary_limit = anti_primary_coordinate_picker(&template_switch.downstream_limit);
        let anti_primary_extended_limit = anti_primary_limit.max(
            anti_primary
                .chars()
                .count()
                .min(template_switch.sp2_secondary_offset + STREAM_PADDING),
        );

        let primary_sp1 = primary_coordinate_picker(&template_switch.sp1_offset);
        let anti_primary_sp1 = anti_primary_coordinate_picker(&template_switch.sp1_offset);
        let primary_sp4 = primary_coordinate_picker(&template_switch.sp4_offset);
        let anti_primary_sp4 = anti_primary_coordinate_picker(&template_switch.sp4_offset);

        debug!("Primary offset: {primary_offset}");
        debug!("SP1 primary offset: {primary_sp1}");
        debug!("SP4 primary offset: {primary_sp4}");
        debug!("Primary limit: {primary_limit}");
        debug!("Primary len: {}", primary.len());
        debug!("Anti-primary extended offset: {anti_primary_extended_offset}");
        debug!("Anti-primary offset: {anti_primary_offset}");
        debug!("SP1 anti-primary offset: {anti_primary_sp1}");
        debug!("SP4 anti-primary offset: {anti_primary_sp4}");
        debug!("Anti-primary limit: {anti_primary_limit}");
        debug!("Anti-primary extended limit: {anti_primary_extended_limit}");
        debug!("Anti-primary len: {}", anti_primary.len());
        debug!(
            "SP2 primary offset: {}",
            template_switch.sp2_secondary_offset
        );
        debug!(
            "SP3 primary offset: {}",
            template_switch.sp3_secondary_offset
        );

        let ts_reverse: String = primary[primary_coordinate_picker(&template_switch.sp1_offset)
            ..primary_coordinate_picker(&template_switch.sp4_offset)]
            .chars()
            .rev()
            .collect();
        let ts_reverse_alignment: Vec<_> = template_switch
            .template_switch
            .iter()
            .copied()
            .rev()
            .collect();

        debug!("Creating renderer");
        let mut renderer = MultipairAlignmentRenderer::new(
            anti_primary_forward_name.clone(),
            &anti_primary[anti_primary_extended_offset..anti_primary_extended_limit],
        );
        debug!("Adding complement primary");
        renderer.add_aligned_sequence(
            &anti_primary_forward_name,
            0,
            anti_primary_reverse_name.clone(),
            &anti_primary_c[anti_primary_extended_offset..anti_primary_extended_limit],
            &[(
                anti_primary_extended_limit - anti_primary_extended_offset,
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );

        debug!("Adding F1");
        renderer.add_aligned_sequence(
            &anti_primary_forward_name,
            anti_primary_offset - anti_primary_extended_offset,
            f1_name.clone(),
            &primary[primary_offset..primary_coordinate_picker(&template_switch.sp1_offset)],
            &template_switch.upstream,
            true,
            invert_alignment,
        );
        debug!("Adding F3");
        renderer.add_aligned_sequence(
            &anti_primary_forward_name,
            anti_primary_coordinate_picker(&template_switch.sp4_offset)
                - anti_primary_extended_offset,
            f3_name.clone(),
            &primary[primary_coordinate_picker(&template_switch.sp4_offset)..primary_limit],
            &template_switch.downstream,
            true,
            invert_alignment,
        );

        debug!("Adding F2");
        renderer.add_aligned_sequence(
            &anti_primary_reverse_name,
            template_switch.sp3_secondary_offset - anti_primary_extended_offset,
            f2_name.clone(),
            &ts_reverse,
            &ts_reverse_alignment,
            true,
            invert_alignment,
        );

        renderer
            .render(
                stdout(),
                &[
                    f1_name,
                    f3_name,
                    anti_primary_forward_name,
                    anti_primary_reverse_name,
                    f2_name,
                ],
            )
            .unwrap();
    }
}
