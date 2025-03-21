use std::{
    fs::File,
    io::{Read, stdout},
    iter,
    path::PathBuf,
};

use alignment_stream::{AlignmentCoordinates, AlignmentStream};
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
use svg::create_ts_svg;

mod alignment_stream;
mod mutlipair_alignment_renderer;
mod parse_template_switches;
mod svg;

#[derive(Parser)]
pub struct Cli {
    #[clap(long, short = 'l', default_value = "info")]
    log_level: LevelFilter,

    /// Path to a toml output file of `tsalign align`.
    ///
    /// This is the input from which template switches are visualised.
    #[clap(long, short = 'i')]
    input: PathBuf,

    /// Path to a toml output file of `tsalign align` created with the `--no-ts` flag.
    ///
    /// This input can be used to provide an alternative alignment without template switches.
    /// It is expected to contain no template switches, otherwise it will be rejected.
    #[clap(long, short = 'n')]
    no_ts_input: Option<PathBuf>,

    /// Create an SVG image file showing the complete alignment with all template switches.
    #[clap(long, short = 's')]
    svg: Option<PathBuf>,

    /// If set together with --svg, then the SVG image is rendered to a PNG image as well.
    #[clap(long, short = 'p')]
    png: bool,
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
    let mut buffer = String::new();
    File::open(cli.input)
        .unwrap_or_else(|error| panic!("Error opening input file: {error}"))
        .read_to_string(&mut buffer)
        .unwrap_or_else(|error| panic!("Error reading input file: {error}"));
    let result =
        toml::from_str(&buffer).unwrap_or_else(|error| panic!("Error parsing input file: {error}"));

    let no_ts_result = cli.no_ts_input.as_ref().map(|no_ts_input| {
        info!("Reading tsalign no-ts output toml file {no_ts_input:?}");
        buffer.clear();
        File::open(no_ts_input)
            .unwrap_or_else(|error| panic!("Error opening input file: {error}"))
            .read_to_string(&mut buffer)
            .unwrap_or_else(|error| panic!("Error reading input file: {error}"));
        toml::from_str(&buffer).unwrap_or_else(|error| panic!("Error parsing input file: {error}"))
    });

    show_template_switches(&result, &no_ts_result);

    if let Some(svg) = cli.svg.as_ref() {
        create_ts_svg(svg, &result, &no_ts_result, cli.png);
    }
}

fn show_template_switches(
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
) {
    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        warn!("Alignment was aborted early, no template switches present");
        return;
    };

    info!("CIGAR: {} (Cost: {:.0})", result.cigar(), statistics.cost);
    if let Some(no_ts_result) = no_ts_result.as_ref() {
        info!(
            "No-ts CIGAR: {} (Cost: {:.0})",
            no_ts_result.cigar(),
            no_ts_result.statistics().cost
        )
    }

    debug!("Collecting template switches");
    let template_switches = parse_template_switches::parse(alignment);
    info!("Found {} template switches", template_switches.len());

    for (index, template_switch) in template_switches.iter().enumerate() {
        info!("Showing template switch {}", index + 1);
        show_template_switch(template_switch, &statistics.sequences, no_ts_result);
    }
}

fn show_template_switch(
    template_switch: &TSShow<AlignmentType>,
    sequences: &SequencePair,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
) {
    // println!("Showing template switch\n{template_switch:?}");

    let reference = &sequences.reference;
    let reference_c: String = sequences.reference_rc.chars().rev().collect();
    let query = &sequences.query;
    let query_c: String = sequences.query_rc.chars().rev().collect();

    let (
        primary_label,
        primary_name,
        primary,
        primary_c,
        primary_coordinate_picker,
        anti_primary_label,
        anti_primary_name,
        anti_primary,
        anti_primary_c,
        anti_primary_coordinate_picker,
        invert_alignment,
    ) = match template_switch.primary {
        TemplateSwitchPrimary::Reference => (
            "Parent".to_string(),
            &sequences.reference_name,
            reference,
            &reference_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| {
                alignment_coordinates.reference()
            }) as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            "Child".to_string(),
            &sequences.query_name,
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| alignment_coordinates.query())
                as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            true,
        ),
        TemplateSwitchPrimary::Query => (
            "Child".to_string(),
            &sequences.query_name,
            query,
            &query_c,
            Box::new(|alignment_coordinates: &AlignmentCoordinates| alignment_coordinates.query())
                as Box<dyn Fn(&AlignmentCoordinates) -> usize>,
            "Parent".to_string(),
            &sequences.reference_name,
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

    let primary_reverse_label = format!("{primary_label}R");
    let anti_primary_forward_label = format!("{anti_primary_label}F");
    let anti_primary_reverse_label = format!("{anti_primary_label}R");
    let f1_label = format!("{primary_label}1");
    let f2_label = format!("{primary_label}2");
    let f3_label = format!("{primary_label}3");

    let primary_offset = primary_coordinate_picker(&template_switch.upstream_offset);
    let primary_limit = primary_coordinate_picker(&template_switch.downstream_limit);
    let anti_primary_offset = anti_primary_coordinate_picker(&template_switch.upstream_offset);
    let anti_primary_limit = anti_primary_coordinate_picker(&template_switch.downstream_limit);

    let primary_sp1 = primary_coordinate_picker(&template_switch.sp1_offset);
    let anti_primary_sp1 = anti_primary_coordinate_picker(&template_switch.sp1_offset);
    let primary_sp4 = primary_coordinate_picker(&template_switch.sp4_offset);
    let anti_primary_sp4 = anti_primary_coordinate_picker(&template_switch.sp4_offset);

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

    debug!("Primary offset: {primary_offset}");
    debug!("SP1 primary offset: {primary_sp1}");
    debug!("SP4 primary offset: {primary_sp4}");
    debug!("Primary limit: {primary_limit}");
    debug!("Primary len: {}", primary.len());
    debug!("Anti-primary offset: {anti_primary_offset}");
    debug!("SP1 anti-primary offset: {anti_primary_sp1}");
    debug!("SP4 anti-primary offset: {anti_primary_sp4}");
    debug!("Anti-primary limit: {anti_primary_limit}");
    debug!("Anti-primary len: {}", anti_primary.len());
    debug!(
        "SP2 primary offset: {}",
        template_switch.sp2_secondary_offset
    );
    debug!(
        "SP3 primary offset: {}",
        template_switch.sp3_secondary_offset
    );

    if primary_equals_secondary {
        let primary_extended_offset = primary_offset.min(
            template_switch
                .sp3_secondary_offset
                .saturating_sub(STREAM_PADDING),
        );
        let primary_extended_limit = primary_limit.max(
            primary
                .chars()
                .count()
                .min(template_switch.sp2_secondary_offset + STREAM_PADDING),
        );

        debug!("Primary extended offset: {primary_extended_offset}");
        debug!("Primary extended limit: {primary_extended_limit}");

        debug!("Creating outside renderer");
        let mut outside_renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_forward_label.clone(),
            anti_primary[anti_primary_offset..anti_primary_limit].chars(),
        );

        debug!("Adding F1");
        outside_renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            0,
            f1_label.clone(),
            primary[primary_offset..primary_coordinate_picker(&template_switch.sp1_offset)].chars(),
            template_switch.upstream.iter().copied(),
            true,
            invert_alignment,
        );
        debug!("Adding F3");
        outside_renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_coordinate_picker(&template_switch.sp4_offset) - anti_primary_offset,
            f3_label.clone(),
            primary[primary_coordinate_picker(&template_switch.sp4_offset)..primary_limit].chars(),
            template_switch.downstream.iter().copied(),
            true,
            invert_alignment,
        );

        debug!("Creating inside renderer");
        let mut inside_renderer = MultipairAlignmentRenderer::new_without_data(
            primary_reverse_label.clone(),
            primary_c[primary_extended_offset..primary_extended_limit].chars(),
        );
        debug!("Adding F2");
        inside_renderer.add_aligned_sequence_without_data(
            &primary_reverse_label,
            template_switch.sp3_secondary_offset - primary_extended_offset,
            f2_label.clone(),
            ts_reverse.chars(),
            ts_reverse_alignment.iter().copied(),
            true,
            false,
        );

        println!("{anti_primary_label}: {anti_primary_name}");
        println!("{primary_label}: {primary_name}");
        println!();
        println!("Switch process:");

        debug!("Rendering");
        outside_renderer
            .render(
                stdout(),
                [&f1_label, &f3_label, &anti_primary_forward_label],
            )
            .unwrap();
        println!();
        inside_renderer
            .render(stdout(), [&primary_reverse_label, &f2_label])
            .unwrap();
    } else {
        let anti_primary_extended_offset = anti_primary_offset.min(
            template_switch
                .sp3_secondary_offset
                .saturating_sub(STREAM_PADDING),
        );
        let anti_primary_extended_limit = anti_primary_limit.max(
            anti_primary
                .chars()
                .count()
                .min(template_switch.sp2_secondary_offset + STREAM_PADDING),
        );

        debug!("Anti-primary extended offset: {anti_primary_extended_offset}");
        debug!("Anti-primary extended limit: {anti_primary_extended_limit}");

        debug!("Creating renderer");
        let mut renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_forward_label.clone(),
            anti_primary[anti_primary_extended_offset..anti_primary_extended_limit].chars(),
        );
        debug!("Adding complement primary");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            0,
            anti_primary_reverse_label.clone(),
            anti_primary_c[anti_primary_extended_offset..anti_primary_extended_limit].chars(),
            [(
                anti_primary_extended_limit - anti_primary_extended_offset,
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );

        debug!("Adding F1");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_offset - anti_primary_extended_offset,
            f1_label.clone(),
            primary[primary_offset..primary_coordinate_picker(&template_switch.sp1_offset)].chars(),
            template_switch.upstream.iter().copied(),
            true,
            invert_alignment,
        );
        debug!("Adding F3");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_forward_label,
            anti_primary_coordinate_picker(&template_switch.sp4_offset)
                - anti_primary_extended_offset,
            f3_label.clone(),
            primary[primary_coordinate_picker(&template_switch.sp4_offset)..primary_limit].chars(),
            template_switch.downstream.iter().copied(),
            true,
            invert_alignment,
        );

        debug!("Adding F2");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_reverse_label,
            template_switch.sp3_secondary_offset - anti_primary_extended_offset,
            f2_label.clone(),
            ts_reverse.chars(),
            ts_reverse_alignment.iter().copied(),
            true,
            false,
        );

        println!("{anti_primary_label}: {anti_primary_name}");
        println!("{primary_label}: {primary_name}");
        println!();
        println!("Switch process:");

        debug!("Rendering");
        renderer
            .render(
                stdout(),
                [
                    &f1_label,
                    &f3_label,
                    &anti_primary_forward_label,
                    &anti_primary_reverse_label,
                    &f2_label,
                ],
            )
            .unwrap();
    }

    if let Some(no_ts_result) = no_ts_result {
        let AlignmentResult::WithTarget {
            alignment: no_ts_alignment,
            ..
        } = no_ts_result
        else {
            warn!("No-ts alignment was aborted early, unable to render");
            return;
        };

        assert!(
            no_ts_alignment.iter().all(|(_, alignment_type)| !matches!(
                alignment_type,
                AlignmentType::TemplateSwitchEntrance { .. }
            )),
            "No-ts alignment must not contain template switches."
        );

        // Find subsequence of no-ts alignment that matches ts alignment interval.
        let mut stream = AlignmentStream::new();
        for alignment_type in no_ts_alignment
            .iter()
            .copied()
            .flat_map(|(multiplicity, alignment_type)| iter::repeat_n(alignment_type, multiplicity))
        {
            if anti_primary_coordinate_picker(&stream.head_coordinates()) >= anti_primary_limit {
                break;
            } else {
                stream.push(1, alignment_type);
            }
        }
        assert_eq!(
            anti_primary_coordinate_picker(&stream.head_coordinates()),
            anti_primary_limit
        );

        while anti_primary_coordinate_picker(&stream.tail_coordinates()) < anti_primary_offset {
            stream.pop_one();
        }
        assert_eq!(
            anti_primary_coordinate_picker(&stream.tail_coordinates()),
            anti_primary_offset
        );

        debug!("Creating no-ts renderer");
        let mut renderer = MultipairAlignmentRenderer::new_without_data(
            anti_primary_label.clone(),
            anti_primary[anti_primary_offset..anti_primary_limit].chars(),
        );

        debug!("Adding primary");
        renderer.add_aligned_sequence_without_data(
            &anti_primary_label,
            0,
            primary_label.clone(),
            primary[primary_coordinate_picker(&stream.tail_coordinates())
                ..primary_coordinate_picker(&stream.head_coordinates())]
                .chars(),
            stream.stream_iter(),
            true,
            invert_alignment,
        );

        println!();
        println!("No-ts alignment:");

        debug!("Rendering");
        renderer
            .render(stdout(), [&anti_primary_label, &primary_label])
            .unwrap();
    } else {
        debug!("No no-ts alignment given, skipping");
    }
}
