use std::{fs::File, io::Read, iter, path::PathBuf};

use clap::Parser;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{a_star_sequences::SequencePair, AlignmentResult},
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
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

fn show_template_switch(template_switch: &TSShow<AlignmentType>, sequences: &SequencePair) {
    // println!("Showing template switch\n{template_switch:?}");

    let (primary_name, primary, primary_rc, primary_upstream_offset) = match template_switch.primary
    {
        TemplateSwitchPrimary::Reference => (
            "Ref".to_string(),
            &sequences.reference,
            &sequences.reference_rc,
            template_switch.upstream_offset.reference(),
        ),
        TemplateSwitchPrimary::Query => (
            "Qry".to_string(),
            &sequences.query,
            &sequences.query_rc,
            template_switch.upstream_offset.query(),
        ),
    };
    let (secondary_name, secondary, secondary_rc, secondary_upstream_offset) =
        match template_switch.secondary {
            TemplateSwitchSecondary::Reference => (
                "Ref".to_string(),
                &sequences.reference,
                &sequences.reference_rc,
                template_switch.upstream_offset.reference(),
            ),
            TemplateSwitchSecondary::Query => (
                "Qry".to_string(),
                &sequences.query,
                &sequences.query_rc,
                template_switch.upstream_offset.query(),
            ),
        };
    let primary_equals_secondary = (template_switch.primary == TemplateSwitchPrimary::Reference
        && template_switch.secondary == TemplateSwitchSecondary::Reference)
        || (template_switch.primary == TemplateSwitchPrimary::Query
            && template_switch.secondary == TemplateSwitchSecondary::Query);

    if primary_equals_secondary {
        todo!()
    } else {
        let mut rendered_primary = String::new();
        let mut rendered_primary_c = String::new();
        let mut rendered_secondary_upstream = String::new();
        let mut rendered_secondary_ts = String::new();
        let mut secondary_ts_spacer = String::new();
        let mut rendered_secondary_downstream = String::new();
        let mut secondary_downstream_spacer = String::new();

        let (primary_upstream_count, secondary_upstream_count) = render_primary_alignment(
            &template_switch.upstream,
            primary.chars().skip(primary_upstream_offset),
            secondary.chars().skip(secondary_upstream_offset),
            &mut rendered_primary,
            &mut rendered_secondary_upstream,
        );

        let primary_len = primary.chars().count();
        render_ts_alignment(
            &template_switch.template_switch,
            primary
                .chars()
                .skip(primary_upstream_offset + primary_upstream_count),
            primary_rc
                .chars()
                .skip(primary_len - template_switch.sp2_offset),
            secondary
                .chars()
                .skip(secondary_upstream_offset + secondary_upstream_count),
            &mut rendered_primary,
            &mut rendered_primary_c,
            &mut rendered_secondary_ts,
        );

        rendered_primary_c = rendered_primary_c.chars().rev().collect();
        rendered_secondary_ts = rendered_secondary_ts.chars().rev().collect();

        println!("{secondary_name}1: L {rendered_secondary_upstream} 1");
        println!(
            "{secondary_name}3:{secondary_downstream_spacer} 4 {rendered_secondary_downstream} R"
        );
        println!("{primary_name}F:   {rendered_primary}");
        println!("{primary_name}R:   {rendered_primary_c}");
        println!("{secondary_name}2:{secondary_ts_spacer} 3 {rendered_secondary_ts} 2");
    }
}

fn render_primary_alignment(
    alignment: &[(usize, AlignmentType)],
    mut reference: impl Iterator<Item = char>,
    mut query: impl Iterator<Item = char>,
    reference_out: &mut impl Extend<char>,
    query_out: &mut impl Extend<char>,
) -> (usize, usize) {
    let mut reference_count = 0;
    let mut query_count = 0;

    for alignment_type in alignment
        .iter()
        .copied()
        .flat_map(|(multiplicity, alignment_type)| iter::repeat(alignment_type).take(multiplicity))
    {
        match alignment_type {
            AlignmentType::PrimaryInsertion | AlignmentType::PrimaryFlankInsertion => {
                reference_out.extend(Some('-'));
                query_out.extend(Some(query.next().unwrap()));
                query_count += 1;
            }
            AlignmentType::PrimaryDeletion | AlignmentType::PrimaryFlankDeletion => {
                reference_out.extend(Some(reference.next().unwrap()));
                query_out.extend(Some('-'));
                reference_count += 1;
            }
            AlignmentType::PrimarySubstitution | AlignmentType::PrimaryFlankSubstitution => {
                reference_out.extend(Some(reference.next().unwrap().to_ascii_lowercase()));
                query_out.extend(Some(query.next().unwrap().to_ascii_lowercase()));
                reference_count += 1;
                query_count += 1;
            }
            AlignmentType::PrimaryMatch | AlignmentType::PrimaryFlankMatch => {
                reference_out.extend(Some(reference.next().unwrap()));
                query_out.extend(Some(query.next().unwrap()));
                reference_count += 1;
                query_count += 1;
            }
            AlignmentType::Root | AlignmentType::PrimaryReentry => { /* Ignore */ }
            AlignmentType::SecondaryInsertion
            | AlignmentType::SecondaryDeletion
            | AlignmentType::SecondarySubstitution
            | AlignmentType::SecondaryMatch
            | AlignmentType::TemplateSwitchEntrance { .. }
            | AlignmentType::TemplateSwitchExit { .. }
            | AlignmentType::SecondaryRoot
            | AlignmentType::PrimaryShortcut { .. } => {
                panic!("Not allowed in primary alignment: {alignment_type:?}")
            }
        }
    }

    (reference_count, query_count)
}

#[expect(clippy::too_many_arguments)]
fn render_ts_alignment(
    alignment: &[(usize, AlignmentType)],
    mut reference: impl Iterator<Item = char>,
    mut reference_rc: impl Iterator<Item = char>,
    mut query: impl Iterator<Item = char>,
    reference_out: &mut impl Extend<char>,
    reference_c_rev_out: &mut impl Extend<char>,
    query_rev_out: &mut impl Extend<char>,
) -> (usize, usize) {
    let mut reference_count = 0;
    let mut query_count = 0;

    for alignment_type in alignment
        .iter()
        .copied()
        .flat_map(|(multiplicity, alignment_type)| iter::repeat(alignment_type).take(multiplicity))
    {
        match alignment_type {
            AlignmentType::SecondaryInsertion => {
                reference_out.extend(Some('-'));
                reference_c_rev_out.extend(Some('-'));
                query_rev_out.extend(Some(query.next().unwrap()));
                query_count += 1;
            }
            AlignmentType::SecondaryDeletion => {
                reference_out.extend(Some(reference.next().unwrap()));
                reference_c_rev_out.extend(Some(reference_rc.next().unwrap()));
                query_rev_out.extend(Some('-'));
                reference_count += 1;
            }
            AlignmentType::SecondarySubstitution => {
                reference_out.extend(Some(reference.next().unwrap()));
                reference_c_rev_out.extend(Some(reference_rc.next().unwrap().to_ascii_lowercase()));
                query_rev_out.extend(Some(query.next().unwrap().to_ascii_lowercase()));
                reference_count += 1;
                query_count += 1;
            }
            AlignmentType::SecondaryMatch => {
                reference_out.extend(Some(reference.next().unwrap()));
                reference_c_rev_out.extend(Some(reference_rc.next().unwrap()));
                query_rev_out.extend(Some(query.next().unwrap()));
                reference_count += 1;
                query_count += 1;
            }
            AlignmentType::Root
            | AlignmentType::PrimaryReentry
            | AlignmentType::TemplateSwitchEntrance { .. }
            | AlignmentType::TemplateSwitchExit { .. }
            | AlignmentType::SecondaryRoot => { /* Ignore */ }
            AlignmentType::PrimaryInsertion
            | AlignmentType::PrimaryFlankInsertion
            | AlignmentType::PrimaryDeletion
            | AlignmentType::PrimaryFlankDeletion
            | AlignmentType::PrimarySubstitution
            | AlignmentType::PrimaryFlankSubstitution
            | AlignmentType::PrimaryMatch
            | AlignmentType::PrimaryFlankMatch
            | AlignmentType::PrimaryShortcut { .. } => {
                panic!("Not allowed in template switch alignment: {alignment_type:?}")
            }
        }
    }

    (reference_count, query_count)
}
