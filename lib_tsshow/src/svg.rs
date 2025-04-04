use std::{collections::BTreeMap, io::Write, iter};

use arrows::{Arrow, ArrowEndpointDirection, add_arrow_defs};
use font::{CharacterData, svg_string, typewriter};
use indexed_str::IndexedStr;
use labelled_sequence::LabelledSequence;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{AlignmentResult, alignment::Alignment},
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info, trace, warn};
use numbers::{Number, NumberAlignment};
use offset_shift::{CharacterIsCopy, OffsetShift};
use svg::{
    Document,
    node::element::{Circle, Group, Text},
};

use crate::{
    error::{Error, Result},
    plain_text::{
        alignment_stream::{AlignmentCoordinates, AlignmentStream},
        mutlipair_alignment_renderer::{Character, MultipairAlignmentRenderer},
    },
    ts_arrangement::{
        TsArrangement, complement::ComplementChar, inner::InnerChar, source::SourceChar,
    },
};

mod arrows;
mod font;
mod indexed_str;
pub mod labelled_sequence;
mod numbers;
mod offset_shift;

const COPY_COLORS: &[&str] = &["#003300", "#006600", "#009900", "#00CC00"];
const COMPLEMENT_SOURCE_HIDDEN_COLOR: &str = "grey";

struct SvgLocation {
    pub x: f32,
    pub y: f32,
}

struct TemplateSwitch {
    label: Option<String>,
    tail: AlignmentCoordinates,
    head: AlignmentCoordinates,
    reference_non_blank_offset: usize,
    reference_non_blank_limit: usize,
    query_non_blank_offset: usize,
    query_non_blank_limit: usize,
    inner_non_blank_length: Option<usize>,
    inner_copy_depths: Vec<CharacterIsCopy>,
    alignment: Alignment<AlignmentType>,
}

impl SvgLocation {
    pub fn as_transform(&self) -> String {
        format!("translate({} {})", self.x, self.y)
    }
}

pub fn create_ts_svg(
    output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    render_arrows: bool,
) -> Result<()> {
    info!("Creating template switch SVG");

    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        return Err(Error::AlignmentHasNoTarget);
    };
    debug!("Alignment: {alignment:?}");

    let reference = &statistics.sequences.reference;
    let query = &statistics.sequences.query;
    let reference_c: String = statistics.sequences.reference_rc.chars().rev().collect();
    let query_c: String = statistics.sequences.query_rc.chars().rev().collect();

    debug!("Computing TS arrangement");
    let mut template_switches = Vec::new();
    let mut ts_arrangement = TsArrangement::new(
        reference.len(),
        query.len(),
        alignment.iter_flat_cloned(),
        &mut template_switches,
    );
    ts_arrangement.remove_empty_columns();
    let ts_arrangement = ts_arrangement;

    let mut ts_group =
        Group::new().set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform());

    let reference = IndexedStr::new(reference);
    let query = IndexedStr::new(query);
    let reference_c = IndexedStr::new(&reference_c);
    let query_c = IndexedStr::new(&query_c);

    let mut y = 0.0;
    for reference_inner in ts_arrangement.reference_inners().iter().rev() {
        ts_group = ts_group.add(svg_string(
            reference_inner.1.iter_values().map(render_inner_char(
                match reference_inner.0.primary {
                    TemplateSwitchPrimary::Reference => &reference,
                    TemplateSwitchPrimary::Query => &query,
                },
            )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        y += typewriter::FONT.character_height;
    }

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference_complement()
            .iter_values()
            .map(render_complement_char(&reference_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference()
            .iter_values()
            .map(render_source_char(&reference)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query()
            .iter_values()
            .map(render_source_char(&query)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query_complement()
            .iter_values()
            .map(render_complement_char(&query_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    y += typewriter::FONT.character_height;

    for query_inner in ts_arrangement.query_inners() {
        ts_group = ts_group.add(svg_string(
            query_inner
                .1
                .iter_values()
                .map(render_inner_char(match query_inner.0.primary {
                    TemplateSwitchPrimary::Reference => &reference,
                    TemplateSwitchPrimary::Query => &query,
                })),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        y += typewriter::FONT.character_height;
    }

    let view_box_width = 20.0 + ts_arrangement.width() as f32 * typewriter::FONT.character_width;
    let view_box_height = 20.0 + y;

    let mut svg = Document::new()
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(ts_group);

    svg = svg.set("viewBox", (0, 0, view_box_width, view_box_height));
    svg::write(output, &svg)?;

    Ok(())
}

fn render_source_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&SourceChar) -> Character<CharacterData> {
    |source_char| match source_char {
        SourceChar::Source {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Hidden { .. } | SourceChar::Blank => Character::new_char_with_default(' '),
    }
}

fn render_complement_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&ComplementChar) -> Character<CharacterData> {
    |complement_char| match complement_char {
        ComplementChar::Complement {
            column,
            lower_case,
            source_hidden,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(
                c,
                if *source_hidden {
                    CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
                } else {
                    Default::default()
                },
            )
        }
        ComplementChar::Gap { source_hidden } => Character::new_char(
            '-',
            if *source_hidden {
                CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
            } else {
                Default::default()
            },
        ),
        ComplementChar::Hidden { .. } | ComplementChar::Blank => {
            Character::new_char(' ', Default::default())
        }
    }
}

fn render_inner_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&InnerChar) -> Character<CharacterData> {
    |inner_char| match inner_char {
        InnerChar::Inner {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };
            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Blank => Character::new_char(' ', Default::default()),
    }
}

fn copy_color(copy_depth: &Option<usize>) -> impl ToString {
    if let Some(copy_depth) = copy_depth {
        COPY_COLORS[copy_depth % COPY_COLORS.len()]
    } else {
        "black"
    }
}

pub fn create_ts_svg_old(
    output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    render_arrows: bool,
) -> Result<()> {
    info!("Creating template switch SVG");

    debug!("Rendering ts alignment");
    let AlignmentResult::WithTarget { alignment, .. } = result else {
        return Err(Error::AlignmentHasNoTarget);
    };
    debug!("Alignment: {alignment:?}");

    let mut has_secondary_reference_ts = false;
    let mut has_secondary_query_ts = false;
    let mut has_self_reference_ts = false;
    let mut has_self_query_ts = false;
    let mut has_negative_anti_primary_gap = false;

    for (primary, secondary) in
        alignment
            .iter_compact()
            .filter_map(|(_, alignment_type)| match alignment_type {
                AlignmentType::TemplateSwitchEntrance {
                    primary, secondary, ..
                } => Some((primary, secondary)),
                AlignmentType::TemplateSwitchExit { anti_primary_gap } => {
                    if *anti_primary_gap < 0 {
                        has_negative_anti_primary_gap = true;
                    }
                    None
                }
                _ => None,
            })
    {
        match (primary, secondary) {
            (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Reference) => {
                has_self_reference_ts = true;
                has_secondary_reference_ts = true;
            }
            (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Query) => {
                has_secondary_query_ts = true;
            }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Reference) => {
                has_secondary_reference_ts = true;
            }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Query) => {
                has_self_query_ts = true;
                has_secondary_query_ts = true;
            }
        }
    }

    if has_negative_anti_primary_gap {
        warn!("Template switch with negative anti-primary gap found, this may be buggy");
    }

    if !has_secondary_reference_ts
        && !has_secondary_query_ts
        && !has_self_reference_ts
        && !has_self_query_ts
    {
        warn!("No template switches found");
    }

    let reference = &result.statistics().sequences.reference;
    let query = &result.statistics().sequences.query;
    let reference_c: String = result
        .statistics()
        .sequences
        .reference_rc
        .chars()
        .rev()
        .collect();
    let query_c: String = result
        .statistics()
        .sequences
        .query_rc
        .chars()
        .rev()
        .collect();
    let reference = LabelledSequence::new("Reference", reference);
    let query = LabelledSequence::new("Query", query);
    let reference_c = LabelledSequence::new("Reference Complement", &reference_c);
    let query_c = LabelledSequence::new("Query Complement", &query_c);

    debug!("Rendering reference and query");
    let mut raw_arrows = Vec::new();
    let mut stream = AlignmentStream::new();
    let mut offset_stream = AlignmentStream::new();
    let mut template_switches = Vec::new();
    let mut offset_shift = OffsetShift::default();

    let mut renderer = MultipairAlignmentRenderer::new_empty();
    renderer.add_empty_independent_sequence(reference.label_cloned());
    renderer.add_empty_independent_sequence(query.label_cloned());
    if has_secondary_reference_ts {
        renderer.add_empty_independent_sequence(reference_c.label_cloned());
    }
    if has_secondary_query_ts {
        renderer.add_empty_independent_sequence(query_c.label_cloned());
    }

    let mut alignment_iter = alignment.iter_flat().copied().peekable();

    while let Some(alignment_type) = alignment_iter.peek().copied() {
        //trace!("Outer {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            render_inter_ts(
                &mut renderer,
                &mut offset_shift,
                &mut raw_arrows,
                &reference,
                &query,
                has_secondary_reference_ts.then_some(&reference_c),
                has_secondary_query_ts.then_some(&query_c),
                &stream,
            );

            stream.clear();

            // Accumulate TS
            for alignment_type in alignment_iter.by_ref() {
                //trace!("Skipping {alignment_type}");

                stream.push(1, alignment_type);
                offset_stream.push(1, alignment_type);

                if matches!(alignment_type, AlignmentType::TemplateSwitchExit { .. }) {
                    break;
                }
            }

            template_switches.push(render_ts_base(
                &mut renderer,
                &mut offset_shift,
                &mut raw_arrows,
                &reference,
                &query,
                has_secondary_reference_ts.then_some(&reference_c),
                has_secondary_query_ts.then_some(&query_c),
                &stream,
            ));

            stream.clear();
        } else {
            alignment_iter.next().unwrap();
            stream.push(1, alignment_type);
            offset_stream.push(1, alignment_type);
        }
    }

    if !stream.is_empty() {
        render_inter_ts(
            &mut renderer,
            &mut offset_shift,
            &mut raw_arrows,
            &reference,
            &query,
            has_secondary_reference_ts.then_some(&reference_c),
            has_secondary_query_ts.then_some(&query_c),
            &stream,
        );
    }

    debug!("Rendering TSes");
    let mut ts_secondary_r_labels = Vec::new();
    let mut ts_secondary_q_labels = Vec::new();
    for (
        index,
        TemplateSwitch {
            label: ts_label,
            tail,
            head,
            inner_non_blank_length,
            inner_copy_depths,
            alignment,
            ..
        },
    ) in template_switches.iter_mut().enumerate()
    {
        let AlignmentType::TemplateSwitchEntrance {
            primary,
            secondary,
            first_offset,
        } = alignment.iter_flat_cloned().next().unwrap()
        else {
            unreachable!();
        };

        let primary = match primary {
            TemplateSwitchPrimary::Reference => {
                reference.substring(tail.reference(), head.reference())
            }
            TemplateSwitchPrimary::Query => query.substring(tail.query(), head.query()),
        };
        let (secondary_c_label, secondary_tail, label_vec) = match secondary {
            TemplateSwitchSecondary::Reference => (
                reference_c.label(),
                tail.reference(),
                &mut ts_secondary_r_labels,
            ),
            TemplateSwitchSecondary::Query => {
                (query_c.label(), tail.query(), &mut ts_secondary_q_labels)
            }
        };

        let label = format!("{} TS{}", primary.label(), index + 1);
        let inner_sequence = primary.sequence();
        let inner_sequence_r: String = inner_sequence.chars().rev().collect();
        let inner_sequence_length = inner_sequence.chars().count();
        let inner_insertion_count = alignment
            .iter_compact()
            .map(|(multiplicity, alignment_type)| {
                if matches!(alignment_type, AlignmentType::SecondaryInsertion) {
                    multiplicity
                } else {
                    0
                }
            })
            .sum::<usize>();

        trace!("Inner alignment: {}", alignment.reverse().cigar());
        trace!("secondary_tail: {secondary_tail}");
        trace!("inner_insertion_count: {inner_insertion_count}");
        trace!("inner_sequence_length: {inner_sequence_length}");

        let inner_limit = secondary_tail as isize + first_offset;
        let inner_offset: usize = (inner_limit + inner_insertion_count as isize
            - inner_sequence_length as isize)
            .try_into()
            .unwrap();

        renderer.add_aligned_sequence(
            secondary_c_label,
            inner_offset,
            label.clone(),
            inner_sequence_r
                .chars()
                .zip(inner_copy_depths.iter().rev())
                .map(|(c, cic)| {
                    Character::new_char(
                        c,
                        match cic {
                            CharacterIsCopy::Yes { depth } => {
                                CharacterData::new_colored(COPY_COLORS[*depth % COPY_COLORS.len()])
                            }
                            CharacterIsCopy::No => Default::default(),
                        },
                    )
                }),
            Default::default,
            Default::default,
            alignment.iter_flat_cloned().skip(1).rev().skip(1),
            true,
            false,
        );

        *inner_non_blank_length = Some(renderer.sequence(&label).len_without_blanks());
        *ts_label = Some(label.clone());

        label_vec.push(label);
    }

    debug!("Translating TS base arrows");
    let mut arrows: Vec<_> = raw_arrows
        .into_iter()
        .map(|arrow| arrow.translate_char_gap_columns_to_real_columns(&renderer))
        .collect();

    debug!("Creating TS numbers and arrows");
    let mut numbers = Vec::new();
    for (
        index,
        TemplateSwitch {
            label,
            reference_non_blank_offset,
            reference_non_blank_limit,
            query_non_blank_offset,
            query_non_blank_limit,
            inner_non_blank_length,
            alignment,
            ..
        },
    ) in template_switches.iter().enumerate()
    {
        let label = label.as_ref().unwrap();
        let inner_sequence = renderer.sequence(label);
        let inner_non_blank_length = inner_non_blank_length.unwrap();

        trace!("inner_non_blank_length: {inner_non_blank_length}");
        trace!("inner_sequence: {inner_sequence}");

        let AlignmentType::TemplateSwitchEntrance { primary, .. } =
            alignment.iter_flat_cloned().next().unwrap()
        else {
            unreachable!()
        };

        let (primary_label, primary_non_blank_offset, primary_non_blank_limit) = match primary {
            TemplateSwitchPrimary::Reference => (
                reference.label(),
                reference_non_blank_offset,
                reference_non_blank_limit,
            ),
            TemplateSwitchPrimary::Query => {
                (query.label(), query_non_blank_offset, query_non_blank_limit)
            }
        };

        let primary_tail = if *primary_non_blank_offset == 0 {
            0
        } else {
            renderer
                .sequence(primary_label)
                .translate_offset_without_blanks(*primary_non_blank_offset - 1)
                .map(|offset| offset + 1)
                .unwrap_or_else(|| renderer.sequence(primary_label).len())
        };
        let primary_head = renderer
            .sequence(primary_label)
            .translate_offset_without_blanks(*primary_non_blank_limit)
            .unwrap_or_else(|| renderer.sequence(primary_label).len());
        let inner_offset = inner_sequence.translate_offset_without_blanks(0);
        let inner_limit = if inner_non_blank_length == 0 {
            None
        } else {
            inner_sequence
                .translate_offset_without_blanks(inner_non_blank_length - 1)
                .map(|limit| limit + 1)
        };

        assert_eq!(inner_offset.is_some(), inner_limit.is_some());
        let inner_offset = inner_offset.unwrap_or(0);
        let inner_limit = inner_limit.unwrap_or(0);

        let letter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().nth(index).unwrap();
        let number_1 = Number::new(
            format!("{letter}1"),
            primary_tail,
            primary_label.clone(),
            NumberAlignment::Left,
            0.5,
        );
        let number_2 = Number::new(
            format!("{letter}2"),
            inner_limit,
            label.clone(),
            NumberAlignment::Left,
            0.5,
        );
        let number_3 = Number::new(
            format!("{letter}3"),
            inner_offset,
            label.clone(),
            NumberAlignment::Right,
            0.5,
        );
        let number_4 = Number::new(
            format!("{letter}4"),
            primary_head,
            primary_label.clone(),
            NumberAlignment::Right,
            0.5,
        );

        if render_arrows {
            // Arrow 1 -> 2
            let arrow = Arrow::new_curved(
                primary_tail,
                number_1.width(),
                primary_label.clone(),
                ArrowEndpointDirection::Forward,
                inner_limit,
                number_2.width(),
                label.clone(),
                ArrowEndpointDirection::Forward,
            );
            debug!("Adding arrow {arrow}");
            arrows.push(arrow);

            // Arrow 3 -> 4
            let arrow = Arrow::new_curved(
                inner_offset,
                number_3.width(),
                label.clone(),
                ArrowEndpointDirection::Backward,
                primary_head,
                number_4.width(),
                primary_label.clone(),
                ArrowEndpointDirection::Backward,
            );
            debug!("Adding arrow {arrow}");
            arrows.push(arrow);
        }

        debug!("Adding number {number_1}");
        numbers.push(number_1);
        debug!("Adding number {number_2}");
        numbers.push(number_2);
        debug!("Adding number {number_3}");
        numbers.push(number_3);
        debug!("Adding number {number_4}");
        numbers.push(number_4);
    }

    debug!("Rendering SVG characters");
    let mut rq_group =
        Group::new().set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform());
    let mut y = 0.0;
    let mut rows = BTreeMap::new();

    if has_secondary_reference_ts {
        for ts_label in &ts_secondary_r_labels {
            rq_group = rq_group.add(svg_string(
                renderer.sequence(ts_label).iter(),
                &SvgLocation { x: 0.0, y },
                &typewriter::FONT,
            ));
            rows.insert(ts_label.clone(), y);
            y += typewriter::FONT.character_height;
        }

        rq_group = rq_group.add(svg_string(
            renderer.sequence(reference_c.label()).iter(),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(reference_c.label_cloned(), y);
        y += typewriter::FONT.character_height;
    }

    rq_group = rq_group.add(svg_string(
        renderer.sequence(reference.label()).iter(),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(reference.label_cloned(), y);
    y += typewriter::FONT.character_height;

    rq_group = rq_group.add(svg_string(
        renderer.sequence(query.label()).iter(),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(query.label_cloned(), y);
    y += typewriter::FONT.character_height;

    if has_secondary_query_ts {
        rq_group = rq_group.add(svg_string(
            renderer.sequence(query_c.label()).iter(),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(query_c.label_cloned(), y);
        y += typewriter::FONT.character_height;

        for ts_label in &ts_secondary_q_labels {
            rq_group = rq_group.add(svg_string(
                renderer.sequence(ts_label).iter(),
                &SvgLocation { x: 0.0, y },
                &typewriter::FONT,
            ));
            rows.insert(ts_label.clone(), y);
            y += typewriter::FONT.character_height;
        }
    }

    trace!("After ts y: {y}");

    debug!("Rendering SVG arrows");
    for arrow in &arrows {
        rq_group = rq_group.add(arrow.render(&rows));
    }

    debug!("Rendering numbers");
    for number in &numbers {
        rq_group = rq_group.add(number.render(&rows));
    }

    let mut view_box_width = offset_stream.len() as f32 * typewriter::FONT.character_width + 20.0;
    let mut view_box_height = y + typewriter::FONT.character_height + 20.0;

    let mut svg = Document::new()
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(rq_group);
    svg = add_arrow_defs(svg);

    if let Some(no_ts_result) = no_ts_result {
        debug!("Rendering no-ts alignment");
        let AlignmentResult::WithTarget {
            alignment: no_ts_alignment,
            ..
        } = no_ts_result
        else {
            warn!("No-ts alignment was aborted early, unable to render");
            return Err(Error::NoTsAlignmentHasNoTarget);
        };

        assert!(
            no_ts_alignment
                .iter_compact()
                .all(|(_, alignment_type)| !matches!(
                    alignment_type,
                    AlignmentType::TemplateSwitchEntrance { .. }
                )),
            "No-ts alignment must not contain template switches."
        );

        debug!("Creating no-ts renderer");
        let mut renderer = MultipairAlignmentRenderer::new_without_data(
            reference.label_cloned(),
            no_ts_result.statistics().sequences.reference.chars(),
        );

        debug!("Adding query");
        renderer.add_aligned_sequence_without_data(
            reference.label(),
            0,
            query.label_cloned(),
            no_ts_result.statistics().sequences.query.chars(),
            no_ts_alignment.iter_flat_cloned(),
            true,
            false,
        );

        debug!("Rendering SVG");
        let reference_svg = svg_string(
            renderer.sequence(reference.label()).iter(),
            &SvgLocation { x: 0.0, y: 0.0 },
            &typewriter::FONT,
        );
        let query_svg = svg_string(
            renderer.sequence(query.label()).iter(),
            &SvgLocation {
                x: 0.0,
                y: 1.0 * typewriter::FONT.character_height,
            },
            &typewriter::FONT,
        );

        let group = Group::new()
            .set(
                "transform",
                SvgLocation {
                    x: 10.0,
                    y: view_box_height - 10.0,
                }
                .as_transform(),
            )
            .add(reference_svg)
            .add(query_svg);
        svg = svg.add(group);
        view_box_width = view_box_width.max(
            renderer
                .sequence(reference.label())
                .len()
                .max(renderer.sequence(query.label()).len()) as f32
                * typewriter::FONT.character_width
                + 20.0,
        );
        view_box_height += 3.0 * typewriter::FONT.character_height;
    } else {
        debug!("No no-ts alignment given, skipping");
    }

    svg = svg.set("viewBox", (0, 0, view_box_width, view_box_height));
    svg::write(output, &svg)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn render_inter_ts(
    renderer: &mut MultipairAlignmentRenderer<String, CharacterData>,
    offset_shift: &mut OffsetShift,
    _arrows: &mut impl Extend<Arrow>,
    reference: &LabelledSequence,
    query: &LabelledSequence,
    reference_c: Option<&LabelledSequence>,
    query_c: Option<&LabelledSequence>,
    stream: &AlignmentStream,
) {
    debug!(
        "Rendering inter-ts from RQ {}/{} to RQ {}/{}",
        stream.tail_coordinates().reference(),
        stream.tail_coordinates().query(),
        stream.head_coordinates().reference(),
        stream.head_coordinates().query()
    );
    debug!("Offset shift: {offset_shift}");

    let reference = reference.substring(
        stream.tail_coordinates().reference(),
        stream.head_coordinates().reference(),
    );
    let query = query.substring(
        stream.tail_coordinates().query(),
        stream.head_coordinates().query(),
    );
    let reference_c = reference_c.map(|reference_c| {
        reference_c.substring(
            stream.tail_coordinates().reference(),
            stream.head_coordinates().reference(),
        )
    });
    let query_c = query_c.map(|query_c| {
        query_c.substring(
            stream.tail_coordinates().query(),
            stream.head_coordinates().query(),
        )
    });

    assert!(
        !reference.is_negative_substring(),
        "inter-ts does not do any jumps"
    );
    assert!(
        !query.is_negative_substring(),
        "inter-ts does not do any jumps"
    );
    debug!("Reference length: {}", reference.sequence().chars().count());
    debug!("Query length: {}", query.sequence().chars().count());

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label())
                    .into_iter()
                    .chain([reference.label(), query.label()])
                    .chain(query_c.as_ref().map(|query_c| query_c.label())),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    let mut reference_copy_char_count = 0;
    renderer.extend_sequence(
        reference.label(),
        reference.sequence().chars().map(|c| {
            Character::new_char(
                c,
                match offset_shift.next_reference() {
                    CharacterIsCopy::Yes { depth } => {
                        reference_copy_char_count += 1;
                        CharacterData::new_colored(COPY_COLORS[depth % COPY_COLORS.len()])
                    }
                    CharacterIsCopy::No => Default::default(),
                },
            )
        }),
        Default::default,
    );
    let reference_copy_char_count = reference_copy_char_count;

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label())
                    .into_iter()
                    .chain([reference.label(), query.label()])
                    .chain(query_c.as_ref().map(|query_c| query_c.label())),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    let mut query_copy_char_count = 0;
    renderer.extend_sequence_with_alignment(
        reference.label(),
        query.label(),
        (stream.tail_coordinates().reference() as isize + offset_shift.reference)
            .try_into()
            .unwrap(),
        query.sequence().chars().map(|c| {
            Character::new_char(
                c,
                match offset_shift.next_query() {
                    CharacterIsCopy::Yes { depth } => {
                        query_copy_char_count += 1;
                        CharacterData::new_colored(COPY_COLORS[depth % COPY_COLORS.len()])
                    }
                    CharacterIsCopy::No => Default::default(),
                },
            )
        }),
        Default::default,
        Default::default,
        stream.stream_iter_flat(),
        true,
        false,
    );
    let query_copy_char_count = query_copy_char_count;

    if let Some(reference_c) = reference_c.as_ref() {
        if reference_copy_char_count < reference_c.sequence().chars().count() {
            trace!("Current renderer content:\n{}", {
                let mut out = Vec::new();
                renderer
                    .render_without_names(
                        &mut out,
                        [reference_c.label(), reference.label(), query.label()]
                            .into_iter()
                            .chain(query_c.as_ref().map(|query_c| query_c.label())),
                    )
                    .unwrap();
                String::from_utf8(out).unwrap()
            });

            renderer.extend_sequence_with_alignment_and_default_data(
                reference.label(),
                reference_c.label(),
                usize::try_from(
                    stream.tail_coordinates().reference() as isize + offset_shift.reference,
                )
                .unwrap()
                    + reference_copy_char_count,
                reference_c
                    .sequence()
                    .chars()
                    .skip(reference_copy_char_count),
                iter::repeat_n(
                    AlignmentType::PrimaryMatch,
                    reference_c.sequence().chars().count(),
                )
                .skip(reference_copy_char_count),
                false,
                false,
            );
        }
    }

    if let Some(query_c) = query_c {
        if query_copy_char_count < query_c.sequence().chars().count() {
            trace!("Current renderer content:\n{}", {
                let mut out = Vec::new();
                renderer
                    .render_without_names(
                        &mut out,
                        reference_c
                            .as_ref()
                            .map(|reference_c| reference_c.label())
                            .into_iter()
                            .chain([reference.label(), query.label(), query_c.label()]),
                    )
                    .unwrap();
                String::from_utf8(out).unwrap()
            });

            renderer.extend_sequence_with_alignment_and_default_data(
                query.label(),
                query_c.label(),
                usize::try_from(stream.tail_coordinates().query() as isize + offset_shift.query)
                    .unwrap()
                    + query_copy_char_count,
                query_c.sequence().chars().skip(query_copy_char_count),
                iter::repeat_n(
                    AlignmentType::PrimaryMatch,
                    query_c.sequence().chars().count(),
                )
                .skip(query_copy_char_count),
                false,
                false,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn render_ts_base(
    renderer: &mut MultipairAlignmentRenderer<String, CharacterData>,
    offset_shift: &mut OffsetShift,
    arrows: &mut impl Extend<Arrow>,
    reference: &LabelledSequence,
    query: &LabelledSequence,
    reference_c: Option<&LabelledSequence>,
    query_c: Option<&LabelledSequence>,
    stream: &AlignmentStream,
) -> TemplateSwitch {
    debug!(
        "Rendering ts base {} from RQ {}/{} to RQ {}/{}",
        stream.stream_iter().next().unwrap().1,
        stream.tail_coordinates().reference(),
        stream.tail_coordinates().query(),
        stream.head_coordinates().reference(),
        stream.head_coordinates().query()
    );
    debug!("Offset shift: {offset_shift}");

    let reference_non_blank_offset = renderer.sequence(reference.label()).len_without_blanks();
    let query_non_blank_offset = renderer.sequence(reference.label()).len_without_blanks();

    let reference = reference.substring(
        stream.tail_coordinates().reference(),
        stream.head_coordinates().reference(),
    );
    let query = query.substring(
        stream.tail_coordinates().query(),
        stream.head_coordinates().query(),
    );
    let reference_c = reference_c.map(|reference_c| {
        reference_c.substring(
            stream.tail_coordinates().reference(),
            stream.head_coordinates().reference(),
        )
    });
    let query_c = query_c.map(|query_c| {
        query_c.substring(
            stream.tail_coordinates().query(),
            stream.head_coordinates().query(),
        )
    });

    let mut alignment = stream.stream_iter();
    let (
        _,
        AlignmentType::TemplateSwitchEntrance {
            primary: ts_primary,
            ..
        },
    ) = alignment.next().unwrap()
    else {
        panic!("Wrong template switch entrance");
    };
    let (_, AlignmentType::TemplateSwitchExit { anti_primary_gap }) = alignment.last().unwrap()
    else {
        panic!("Wrong template switch exit")
    };

    let (primary, anti_primary, primary_c, anti_primary_c, primary_offset, anti_primary_offset) =
        match ts_primary {
            TemplateSwitchPrimary::Reference => (
                &reference,
                &query,
                &reference_c,
                &query_c,
                usize::try_from(
                    stream.tail_coordinates().reference() as isize + offset_shift.reference,
                )
                .unwrap(),
                (stream.tail_coordinates().query() as isize + offset_shift.query)
                    .try_into()
                    .unwrap(),
            ),
            TemplateSwitchPrimary::Query => (
                &query,
                &reference,
                &query_c,
                &reference_c,
                (stream.tail_coordinates().query() as isize + offset_shift.query)
                    .try_into()
                    .unwrap(),
                (stream.tail_coordinates().reference() as isize + offset_shift.reference)
                    .try_into()
                    .unwrap(),
            ),
        };

    assert!(!primary.is_negative_substring());
    let primary_copy_depths: Vec<_> = primary
        .sequence()
        .chars()
        .map(|_| offset_shift.next_primary(ts_primary))
        .collect();
    let primary_copied_char_count = primary_copy_depths
        .iter()
        .filter(|cic| matches!(cic, CharacterIsCopy::Yes { .. }))
        .count();

    if !anti_primary.is_negative_substring() {
        let mut anti_primary_copied_count = 0;
        renderer.extend_sequence(
            anti_primary.label(),
            anti_primary.sequence().chars().map(|c| {
                Character::new_char(c, {
                    match offset_shift.next_anti_primary(ts_primary) {
                        CharacterIsCopy::Yes { depth } => {
                            anti_primary_copied_count += 1;
                            CharacterData::new_colored(COPY_COLORS[depth % COPY_COLORS.len()])
                        }
                        CharacterIsCopy::No => Default::default(),
                    }
                })
            }),
            Default::default,
        );
        let anti_primary_copied_char_count = anti_primary_copied_count;

        trace!("Current renderer content:\n{}", {
            let mut out = Vec::new();
            renderer
                .render_without_names(
                    &mut out,
                    reference_c
                        .as_ref()
                        .map(|reference_c| reference_c.label())
                        .into_iter()
                        .chain([reference.label(), query.label()])
                        .chain(query_c.as_ref().map(|query_c| query_c.label())),
                )
                .unwrap();
            String::from_utf8(out).unwrap()
        });

        renderer.extend_sequence_with_alignment(
            anti_primary.label(),
            primary.label(),
            anti_primary_offset,
            anti_primary.sequence().chars().map(|c| {
                Character::new_char(
                    c,
                    CharacterData {
                        color: "white".to_string(),
                    },
                )
            }),
            Default::default,
            || unreachable!(),
            iter::repeat_n(
                AlignmentType::PrimaryMatch,
                anti_primary.sequence().chars().count(),
            ),
            false,
            false,
        );

        if let Some(anti_primary_c) = anti_primary_c {
            trace!("Current renderer content:\n{}", {
                let mut out = Vec::new();
                renderer
                    .render_without_names(
                        &mut out,
                        reference_c
                            .as_ref()
                            .map(|reference_c| reference_c.label())
                            .into_iter()
                            .chain([reference.label(), query.label()])
                            .chain(query_c.as_ref().map(|query_c| query_c.label())),
                    )
                    .unwrap();
                String::from_utf8(out).unwrap()
            });

            for (offset, character) in anti_primary_c
                .sequence()
                .char_indices()
                .skip(anti_primary_copied_char_count)
            {
                renderer.extend_sequence_with_alignment_and_default_data(
                    anti_primary.label(),
                    anti_primary_c.label(),
                    offset_shift.anti_primary_source_character_to_rendered_character_index(
                        ts_primary,
                        anti_primary_offset + offset,
                    ),
                    iter::once(character),
                    iter::once(AlignmentType::PrimaryMatch),
                    false,
                    false,
                );
            }
        }

        if let Some(primary_c) = &primary_c {
            let extension = primary_c
                .sequence()
                .chars()
                .map(|c| Character::new_char(c, CharacterData::new_colored("grey")));

            let extension_total_length = extension.clone().count();
            let extension_existing_length = anti_primary.sequence().chars().count();
            let extension_additional_length =
                extension_total_length.saturating_sub(extension_existing_length);
            trace!("extension_existing_length: {extension_existing_length}");
            trace!("extension_additional_length: {extension_additional_length}");

            trace!("Current renderer content:\n{}", {
                let mut out = Vec::new();
                renderer
                    .render_without_names(
                        &mut out,
                        reference_c
                            .as_ref()
                            .map(|reference_c| reference_c.label())
                            .into_iter()
                            .chain([reference.label(), query.label()])
                            .chain(query_c.as_ref().map(|query_c| query_c.label())),
                    )
                    .unwrap();
                String::from_utf8(out).unwrap()
            });

            for (offset, character) in extension
                .clone()
                .enumerate()
                .take(extension_existing_length)
                .skip(primary_copied_char_count)
            {
                renderer.extend_sequence_with_alignment(
                    primary.label(),
                    primary_c.label(),
                    offset_shift.primary_source_character_to_rendered_character_index(
                        ts_primary,
                        primary_offset + offset,
                    ),
                    iter::once(character),
                    Default::default,
                    Default::default,
                    iter::once(AlignmentType::PrimaryMatch),
                    false,
                    false,
                );
            }

            for (offset, character) in extension
                .clone()
                .enumerate()
                .skip(extension_existing_length.max(primary_copied_char_count))
            {
                renderer.extend_sequence_with_alignment(
                    primary.label(),
                    primary_c.label(),
                    offset_shift.primary_source_character_to_rendered_character_index(
                        ts_primary,
                        primary_offset + offset,
                    ),
                    iter::once(character),
                    Default::default,
                    Default::default,
                    iter::once(AlignmentType::PrimaryMatch),
                    false,
                    false,
                );
            }

            if extension_additional_length > 0 {
                let column_a = renderer.sequence(anti_primary.label()).len_without_blanks();
                let arrow_a = Arrow::new_skip(column_a, column_a, anti_primary.label().clone());
                debug!("Adding extension arrow {arrow_a}");
                arrows.extend([arrow_a]);

                if let Some(anti_primary_c) = &anti_primary_c {
                    let column_ac = renderer
                        .sequence(anti_primary_c.label())
                        .len_without_blanks();
                    let arrow_ac =
                        Arrow::new_skip(column_ac, column_ac, anti_primary_c.label().clone());
                    debug!("Adding extension arrow {arrow_ac}");
                    arrows.extend([arrow_ac]);
                }
            } else if extension_total_length < extension_existing_length {
                let column_c = renderer.sequence(primary_c.label()).len_without_blanks();
                let arrow_c = Arrow::new_skip(column_c, column_c, primary_c.label().clone());
                debug!("Adding extension arrow {arrow_c}");
                arrows.extend([arrow_c]);
            }
        }
    } else {
        trace!("Anti-primary is negative substring, no need to render anything");
    }

    let length_difference =
        anti_primary_gap - isize::try_from(primary.sequence().chars().count()).unwrap();
    match ts_primary {
        TemplateSwitchPrimary::Reference => {
            offset_shift.reference += length_difference;
            if query.is_negative_substring() {
                offset_shift.query -= isize::try_from(query.sequence().chars().count()).unwrap();
                offset_shift.push_query_copy(query.sequence().chars().count());
            }
        }
        TemplateSwitchPrimary::Query => {
            offset_shift.query += length_difference;
            if reference.is_negative_substring() {
                offset_shift.reference -=
                    isize::try_from(reference.sequence().chars().count()).unwrap();
                offset_shift.push_reference_copy(reference.sequence().chars().count());
            }
        }
    }

    let reference_non_blank_limit = renderer.sequence(reference.label()).len_without_blanks();
    let query_non_blank_limit = renderer.sequence(reference.label()).len_without_blanks();

    TemplateSwitch {
        label: None,
        tail: stream.tail_coordinates(),
        head: stream.head_coordinates(),
        reference_non_blank_offset,
        reference_non_blank_limit,
        query_non_blank_offset,
        query_non_blank_limit,
        inner_non_blank_length: None,
        inner_copy_depths: primary_copy_depths,
        alignment: stream.stream_alignment(),
    }
}

pub fn create_error_svg(output: impl Write, error: Error) -> Result<()> {
    let svg = Document::new()
        .set("viewBox", (0, 0, 1920, 1080))
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(
            Text::new(error.to_string())
                .set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform())
                .set("font", "Sans serif"),
        );

    svg::write(output, &svg)?;
    Ok(())
}
