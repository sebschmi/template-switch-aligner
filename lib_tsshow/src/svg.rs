use core::str;
use std::{collections::BTreeMap, io::Write, iter};

use arrows::{Arrow, ArrowEndpointDirection, add_arrow_defs};
use font::{CharacterData, svg_string, typewriter};
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::{AlignmentResult, alignment::Alignment},
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info, trace, warn};
use numbers::{Number, NumberAlignment};
use svg::{
    Document,
    node::element::{Circle, Group},
};

use crate::plain_text::{
    alignment_stream::{AlignmentCoordinates, AlignmentStream},
    mutlipair_alignment_renderer::{Character, MultipairAlignmentRenderer},
};

mod arrows;
mod font;
mod numbers;

struct SvgLocation {
    pub x: f32,
    pub y: f32,
}

#[derive(Debug, Clone)]
struct LabelledSequence<'a, Label> {
    label: Label,
    sequence: &'a str,
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
    alignment: Alignment<AlignmentType>,
}

#[derive(Debug, Default)]
struct OffsetShift {
    reference: isize,
    query: isize,
}

impl SvgLocation {
    pub fn as_transform(&self) -> String {
        format!("translate({} {})", self.x, self.y)
    }
}

impl<Label: Clone> LabelledSequence<'_, Label> {
    fn substring(&self, offset: usize, limit: usize) -> Self {
        Self {
            label: self.label.clone(),
            sequence: char_substring(self.sequence, offset, limit),
        }
    }
}

pub fn create_ts_svg(
    output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    render_arrows: bool,
) {
    info!("Creating template switch SVG");

    debug!("Rendering ts alignment");
    let AlignmentResult::WithTarget { alignment, .. } = result else {
        warn!("Alignment was aborted early, unable to render");
        return;
    };
    debug!("Alignment: {alignment:?}");

    let mut has_secondary_reference_ts = false;
    let mut has_secondary_query_ts = false;
    let mut has_self_reference_ts = false;
    let mut has_self_query_ts = false;

    for (primary, secondary) in alignment.iter_compact().filter_map(|(_, alignment_type)| {
        if let AlignmentType::TemplateSwitchEntrance {
            primary, secondary, ..
        } = *alignment_type
        {
            Some((primary, secondary))
        } else {
            None
        }
    }) {
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

    if !has_secondary_reference_ts
        && !has_secondary_query_ts
        && !has_self_reference_ts
        && !has_self_query_ts
    {
        warn!("No template switches found");
    }

    let reference_label = "Reference".to_string();
    let query_label = "Query".to_string();
    let reference_c_label = "Reference Complement".to_string();
    let query_c_label = "Query Complement".to_string();

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
    let reference = LabelledSequence {
        label: &reference_label,
        sequence: reference,
    };
    let query = LabelledSequence {
        label: &query_label,
        sequence: query,
    };
    let reference_c = LabelledSequence {
        label: &reference_c_label,
        sequence: &reference_c,
    };
    let query_c = LabelledSequence {
        label: &query_c_label,
        sequence: &query_c,
    };

    debug!("Rendering reference and query");
    let mut arrows = Vec::new();
    let mut stream = AlignmentStream::new();
    let mut offset_stream = AlignmentStream::new();
    let mut template_switches = Vec::new();
    let mut offset_shift = OffsetShift::default();

    let mut renderer = MultipairAlignmentRenderer::new_empty();
    renderer.add_empty_independent_sequence(reference_label.clone());
    renderer.add_empty_independent_sequence(query_label.clone());
    if has_secondary_reference_ts {
        renderer.add_empty_independent_sequence(reference_c_label.clone());
    }
    if has_secondary_query_ts {
        renderer.add_empty_independent_sequence(query_c_label.clone());
    }

    let mut alignment_iter = alignment.iter_flat().copied().peekable();

    while let Some(alignment_type) = alignment_iter.peek().copied() {
        //trace!("Outer {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            render_inter_ts(
                &mut renderer,
                &offset_shift,
                &mut arrows,
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
                &mut arrows,
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
            &offset_shift,
            &mut arrows,
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
                reference_c.label,
                tail.reference(),
                &mut ts_secondary_r_labels,
            ),
            TemplateSwitchSecondary::Query => {
                (query_c.label, tail.query(), &mut ts_secondary_q_labels)
            }
        };

        let label = format!("{} TS{}", primary.label, index + 1);
        let inner_sequence = primary.sequence;
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

        renderer.add_aligned_sequence_with_default_data(
            secondary_c_label,
            inner_offset,
            label.clone(),
            inner_sequence_r.chars(),
            alignment.iter_flat_cloned().skip(1).rev().skip(1),
            true,
            false,
        );

        *inner_non_blank_length = Some(renderer.sequence(&label).len_without_blanks());
        *ts_label = Some(label.clone());

        label_vec.push(label);
    }

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
                &reference_label,
                reference_non_blank_offset,
                reference_non_blank_limit,
            ),
            TemplateSwitchPrimary::Query => {
                (&query_label, query_non_blank_offset, query_non_blank_limit)
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
            renderer.sequence(&reference_c_label).iter(),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(reference_c_label.clone(), y);
        y += typewriter::FONT.character_height;
    }

    rq_group = rq_group.add(svg_string(
        renderer.sequence(&reference_label).iter(),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(reference_label.clone(), y);
    y += typewriter::FONT.character_height;

    rq_group = rq_group.add(svg_string(
        renderer.sequence(&query_label).iter(),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(query_label.clone(), y);
    y += typewriter::FONT.character_height;

    if has_secondary_query_ts {
        rq_group = rq_group.add(svg_string(
            renderer.sequence(&query_c_label).iter(),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(query_c_label.clone(), y);
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
            return;
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
            reference_label.clone(),
            no_ts_result.statistics().sequences.reference.chars(),
        );

        debug!("Adding query");
        renderer.add_aligned_sequence_without_data(
            &reference_label,
            0,
            query_label.clone(),
            no_ts_result.statistics().sequences.query.chars(),
            no_ts_alignment.iter_flat_cloned(),
            true,
            false,
        );

        debug!("Rendering SVG");
        let reference = svg_string(
            renderer.sequence(&reference_label).iter(),
            &SvgLocation { x: 0.0, y: 0.0 },
            &typewriter::FONT,
        );
        let query = svg_string(
            renderer.sequence(&query_label).iter(),
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
            .add(reference)
            .add(query);
        svg = svg.add(group);
        view_box_width = view_box_width.max(
            renderer
                .sequence(&reference_label)
                .len()
                .max(renderer.sequence(&query_label).len()) as f32
                * typewriter::FONT.character_width
                + 20.0,
        );
        view_box_height += 3.0 * typewriter::FONT.character_height;
    } else {
        debug!("No no-ts alignment given, skipping");
    }

    svg = svg.set("viewBox", (0, 0, view_box_width, view_box_height));
    svg::write(output, &svg).unwrap();
}

#[allow(clippy::too_many_arguments)]
fn render_inter_ts(
    renderer: &mut MultipairAlignmentRenderer<String, CharacterData>,
    offset_shift: &OffsetShift,
    _arrows: &mut impl Extend<Arrow>,
    reference: &LabelledSequence<&String>,
    query: &LabelledSequence<&String>,
    reference_c: Option<&LabelledSequence<&String>>,
    query_c: Option<&LabelledSequence<&String>>,
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

    debug!("Reference length: {}", reference.sequence.chars().count());
    debug!("Query length: {}", query.sequence.chars().count());

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label])
                    .chain(query_c.as_ref().map(|query_c| query_c.label)),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    renderer.extend_sequence_with_default_data(reference.label, reference.sequence.chars());

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label])
                    .chain(query_c.as_ref().map(|query_c| query_c.label)),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    renderer.extend_sequence_with_alignment_and_default_data(
        reference.label,
        query.label,
        (stream.tail_coordinates().reference() as isize + offset_shift.reference)
            .try_into()
            .unwrap(),
        query.sequence.chars(),
        stream.stream_iter_flat(),
        true,
        false,
    );

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label])
                    .chain(query_c.as_ref().map(|query_c| query_c.label)),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    if let Some(reference_c) = reference_c.as_ref() {
        renderer.extend_sequence_with_alignment_and_default_data(
            reference.label,
            reference_c.label,
            (stream.tail_coordinates().reference() as isize + offset_shift.reference)
                .try_into()
                .unwrap(),
            reference_c.sequence.chars(),
            iter::repeat_n(
                AlignmentType::PrimaryMatch,
                reference_c.sequence.chars().count(),
            ),
            false,
            false,
        );
    }

    if let Some(query_c) = query_c {
        trace!("Current renderer content:\n{}", {
            let mut out = Vec::new();
            renderer
                .render_without_names(
                    &mut out,
                    reference_c
                        .as_ref()
                        .map(|reference_c| reference_c.label)
                        .into_iter()
                        .chain([reference.label, query.label, query_c.label]),
                )
                .unwrap();
            String::from_utf8(out).unwrap()
        });

        renderer.extend_sequence_with_alignment_and_default_data(
            query.label,
            query_c.label,
            (stream.tail_coordinates().query() as isize + offset_shift.query)
                .try_into()
                .unwrap(),
            query_c.sequence.chars(),
            iter::repeat_n(
                AlignmentType::PrimaryMatch,
                query_c.sequence.chars().count(),
            ),
            false,
            false,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn render_ts_base(
    renderer: &mut MultipairAlignmentRenderer<String, CharacterData>,
    offset_shift: &mut OffsetShift,
    arrows: &mut impl Extend<Arrow>,
    reference: &LabelledSequence<&String>,
    query: &LabelledSequence<&String>,
    reference_c: Option<&LabelledSequence<&String>>,
    query_c: Option<&LabelledSequence<&String>>,
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

    let reference_non_blank_offset = renderer.sequence(reference.label).len_without_blanks();
    let query_non_blank_offset = renderer.sequence(reference.label).len_without_blanks();

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
            secondary: ts_secondary,
            ..
        },
    ) = alignment.next().unwrap()
    else {
        panic!("Wrong template switch entrance");
    };
    let (_, AlignmentType::TemplateSwitchExit { length_difference }) = alignment.last().unwrap()
    else {
        panic!("Wrong template switch exit")
    };

    let (
        secondary,
        anti_secondary,
        secondary_c,
        anti_secondary_c,
        secondary_offset,
        anti_secondary_offset,
    ) = match ts_secondary {
        TemplateSwitchSecondary::Reference => (
            reference.clone(),
            query.clone(),
            reference_c.clone(),
            query_c.clone(),
            (stream.tail_coordinates().reference() as isize + offset_shift.reference)
                .try_into()
                .unwrap(),
            (stream.tail_coordinates().query() as isize + offset_shift.query)
                .try_into()
                .unwrap(),
        ),
        TemplateSwitchSecondary::Query => (
            query.clone(),
            reference.clone(),
            query_c.clone(),
            reference_c.clone(),
            (stream.tail_coordinates().query() as isize + offset_shift.query)
                .try_into()
                .unwrap(),
            (stream.tail_coordinates().reference() as isize + offset_shift.reference)
                .try_into()
                .unwrap(),
        ),
    };

    renderer.extend_sequence_with_default_data(secondary.label, secondary.sequence.chars());

    trace!("Current renderer content:\n{}", {
        let mut out = Vec::new();
        renderer
            .render_without_names(
                &mut out,
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label])
                    .chain(query_c.as_ref().map(|query_c| query_c.label)),
            )
            .unwrap();
        String::from_utf8(out).unwrap()
    });

    renderer.extend_sequence_with_alignment(
        secondary.label,
        anti_secondary.label,
        secondary_offset,
        secondary.sequence.chars().map(|c| {
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
            secondary.sequence.chars().count(),
        ),
        false,
        false,
    );

    if let Some(secondary_c) = &secondary_c {
        renderer.extend_sequence_with_alignment_and_default_data(
            secondary.label,
            secondary_c.label,
            secondary_offset,
            secondary_c.sequence.chars(),
            iter::repeat_n(
                AlignmentType::PrimaryMatch,
                secondary_c.sequence.chars().count(),
            ),
            false,
            false,
        );
    }

    if let Some(anti_secondary_c) = anti_secondary_c {
        let extension = anti_secondary_c.sequence.chars().map(|c| {
            Character::new_char(
                c,
                CharacterData {
                    color: "grey".to_string(),
                },
            )
        });
        let extension_existing_length = secondary.sequence.chars().count();
        let extension_additional_length = extension
            .clone()
            .count()
            .saturating_sub(extension_existing_length);
        trace!("extension_existing_length: {extension_existing_length}");
        trace!("extension_additional_length: {extension_existing_length}");

        trace!("Current renderer content:\n{}", {
            let mut out = Vec::new();
            renderer
                .render_without_names(
                    &mut out,
                    reference_c
                        .as_ref()
                        .map(|reference_c| reference_c.label)
                        .into_iter()
                        .chain([reference.label, query.label])
                        .chain(query_c.as_ref().map(|query_c| query_c.label)),
                )
                .unwrap();
            String::from_utf8(out).unwrap()
        });

        renderer.extend_sequence_with_alignment(
            anti_secondary.label,
            anti_secondary_c.label,
            anti_secondary_offset,
            extension.clone().take(extension_existing_length),
            Default::default,
            Default::default,
            iter::repeat_n(
                AlignmentType::PrimaryMatch,
                extension_existing_length.min(extension.clone().count()),
            ),
            false,
            false,
        );
        let left_extension_column = renderer.column_width();
        renderer.extend_sequence(
            anti_secondary_c.label,
            extension.skip(extension_existing_length),
            Default::default,
        );
        let right_extension_column = renderer.column_width();

        if extension_additional_length > 0 {
            if let Some(secondary_c) = &secondary_c {
                debug!("Adding extension arrow");
                arrows.extend([Arrow::new_skip(
                    left_extension_column,
                    right_extension_column,
                    secondary_c.label.clone(),
                )]);
            }
        }
    }

    match (ts_primary, ts_secondary) {
        (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Reference) => {
            offset_shift.query -= length_difference
        }
        (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Query) => {
            offset_shift.reference += length_difference
        }
        (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Reference) => {
            offset_shift.query += length_difference
        }
        (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Query) => {
            offset_shift.reference -= length_difference
        }
    }

    let reference_non_blank_limit = renderer.sequence(reference.label).len_without_blanks();
    let query_non_blank_limit = renderer.sequence(reference.label).len_without_blanks();

    TemplateSwitch {
        label: None,
        tail: stream.tail_coordinates(),
        head: stream.head_coordinates(),
        reference_non_blank_offset,
        reference_non_blank_limit,
        query_non_blank_offset,
        query_non_blank_limit,
        inner_non_blank_length: None,
        alignment: stream.stream_alignment(),
    }
}

fn char_substring(string: &str, offset: usize, limit: usize) -> &str {
    //trace!("Taking substring {offset}..{limit} of {string}");

    let mut indices = string
        .char_indices()
        .map(|(index, _)| index)
        .chain(Some(string.len()));
    let byte_offset = indices.by_ref().nth(offset).unwrap_or_else(|| {
        panic!(
            "The string contains {} characters, but the offset is {offset}",
            string.chars().count()
        )
    });

    if offset == limit {
        &string[byte_offset..byte_offset]
    } else {
        let byte_limit = indices.nth(limit - offset - 1).unwrap_or_else(|| {
            panic!(
                "The string contains {} characters, but the limit is {limit}",
                string.chars().count()
            )
        });
        &string[byte_offset..byte_limit]
    }
}

impl std::fmt::Display for OffsetShift {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RQ {}/{}", self.reference, self.query)
    }
}
