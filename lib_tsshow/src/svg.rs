use std::{io::stdout, iter};

use font::{svg_string, CharacterData, CHARACTER_HEIGHT, CHARACTER_WIDTH};
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::AlignmentResult,
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info, trace, warn};
use resvg::{
    tiny_skia,
    usvg::{self, Transform},
};
use svg::{
    node::element::{Circle, Group},
    Document,
};

use crate::plain_text::{
    alignment_stream::{AlignmentCoordinates, AlignmentStream},
    mutlipair_alignment_renderer::{Character, MultipairAlignmentRenderer},
};

mod font;

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
    tail: AlignmentCoordinates,
    head: AlignmentCoordinates,
    alignment: Vec<(usize, AlignmentType)>,
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
    output: impl AsRef<std::path::Path>,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    render_png: bool,
) {
    let output = output.as_ref();
    info!("Creating template switch SVG at {output:?}");

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

    for (primary, secondary) in alignment.iter().filter_map(|(_, alignment_type)| {
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

    let mut alignment_iter = alignment
        .iter()
        .copied()
        .flat_map(|(multiplicity, alignment_type)| {
            iter::repeat_n(
                alignment_type,
                if matches!(
                    alignment_type,
                    AlignmentType::TemplateSwitchEntrance { .. }
                        | AlignmentType::TemplateSwitchExit { .. }
                ) {
                    1
                } else {
                    multiplicity
                },
            )
        })
        .peekable();

    while let Some(alignment_type) = alignment_iter.peek().copied() {
        //trace!("Outer {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            render_inter_ts(
                &mut renderer,
                &offset_shift,
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

            template_switches.push(TemplateSwitch {
                tail: stream.tail_coordinates(),
                head: stream.head_coordinates(),
                alignment: stream.stream_vec(),
            });

            render_ts_base(
                &mut renderer,
                &mut offset_shift,
                &reference,
                &query,
                has_secondary_reference_ts.then_some(&reference_c),
                has_secondary_query_ts.then_some(&query_c),
                &stream,
            );

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
            tail,
            head,
            alignment,
        },
    ) in template_switches.iter().enumerate()
    {
        let (
            _,
            AlignmentType::TemplateSwitchEntrance {
                primary,
                secondary,
                first_offset,
            },
        ) = alignment.first().copied().unwrap()
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
            .iter()
            .map(|(multiplicity, alignment_type)| {
                if matches!(alignment_type, AlignmentType::SecondaryInsertion) {
                    *multiplicity
                } else {
                    0
                }
            })
            .sum::<usize>();

        trace!(
            "Inner alignment: {:?}",
            alignment.iter().copied().rev().collect::<Vec<_>>()
        );
        trace!("secondary_tail: {secondary_tail}");
        trace!("inner_insertion_count: {inner_insertion_count}");
        trace!("inner_sequence_length: {inner_sequence_length}");

        let inner_offset: usize = ((secondary_tail + inner_insertion_count) as isize
            + first_offset
            - inner_sequence_length as isize)
            .try_into()
            .unwrap();

        renderer.add_aligned_sequence_with_default_data(
            secondary_c_label,
            inner_offset,
            label.clone(),
            inner_sequence_r.chars(),
            alignment[1..alignment.len() - 1].iter().copied().rev(),
            true,
            false,
        );

        label_vec.push(label);
    }

    debug!("Rendering SVG");
    let mut rq_group =
        Group::new().set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform());
    let mut y = 0.0;

    if has_secondary_reference_ts {
        for ts_label in &ts_secondary_r_labels {
            rq_group = rq_group.add(svg_string(
                renderer.sequence(ts_label).iter(),
                &SvgLocation { x: 0.0, y },
            ));
            y += CHARACTER_HEIGHT;
        }

        rq_group = rq_group.add(svg_string(
            renderer.sequence(&reference_c_label).iter(),
            &SvgLocation { x: 0.0, y },
        ));
        y += CHARACTER_HEIGHT;
    }

    rq_group = rq_group.add(svg_string(
        renderer.sequence(&reference_label).iter(),
        &SvgLocation { x: 0.0, y },
    ));
    y += CHARACTER_HEIGHT;

    rq_group = rq_group.add(svg_string(
        renderer.sequence(&query_label).iter(),
        &SvgLocation { x: 0.0, y },
    ));
    y += CHARACTER_HEIGHT;

    if has_secondary_query_ts {
        rq_group = rq_group.add(svg_string(
            renderer.sequence(&query_c_label).iter(),
            &SvgLocation { x: 0.0, y },
        ));
        y += CHARACTER_HEIGHT;

        for ts_label in &ts_secondary_q_labels {
            rq_group = rq_group.add(svg_string(
                renderer.sequence(ts_label).iter(),
                &SvgLocation { x: 0.0, y },
            ));
            y += CHARACTER_HEIGHT;
        }
    }

    trace!("After ts y: {y}");

    let mut view_box_width = offset_stream.len() as f32 * CHARACTER_WIDTH + 20.0;
    let mut view_box_height = y + CHARACTER_HEIGHT + 20.0;

    let mut svg = Document::new()
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(rq_group);

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
            no_ts_alignment.iter().all(|(_, alignment_type)| !matches!(
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
            no_ts_alignment.iter().copied(),
            true,
            false,
        );

        debug!("Rendering SVG");
        let reference = svg_string(
            renderer.sequence(&reference_label).iter(),
            &SvgLocation { x: 0.0, y: 0.0 },
        );
        let query = svg_string(
            renderer.sequence(&query_label).iter(),
            &SvgLocation {
                x: 0.0,
                y: 1.0 * CHARACTER_HEIGHT,
            },
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
                * CHARACTER_WIDTH
                + 20.0,
        );
        view_box_height += 3.0 * CHARACTER_HEIGHT;
    } else {
        debug!("No no-ts alignment given, skipping");
    }

    svg = svg.set("viewBox", (0, 0, view_box_width, view_box_height));
    svg::save(output, &svg).unwrap();

    if render_png {
        make_png(output);
    }
}

fn render_inter_ts<SequenceName: Eq + Ord>(
    renderer: &mut MultipairAlignmentRenderer<SequenceName, CharacterData>,
    offset_shift: &OffsetShift,
    reference: &LabelledSequence<&SequenceName>,
    query: &LabelledSequence<&SequenceName>,
    reference_c: Option<&LabelledSequence<&SequenceName>>,
    query_c: Option<&LabelledSequence<&SequenceName>>,
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

    renderer
        .render_without_names(
            stdout(),
            reference_c
                .as_ref()
                .map(|reference_c| reference_c.label)
                .into_iter()
                .chain([reference.label, query.label])
                .chain(query_c.as_ref().map(|query_c| query_c.label)),
        )
        .unwrap();

    renderer.extend_sequence_with_default_data(reference.label, reference.sequence.chars());

    renderer
        .render_without_names(
            stdout(),
            reference_c
                .as_ref()
                .map(|reference_c| reference_c.label)
                .into_iter()
                .chain([reference.label, query.label])
                .chain(query_c.as_ref().map(|query_c| query_c.label)),
        )
        .unwrap();

    renderer.extend_sequence_with_alignment_and_default_data(
        reference.label,
        query.label,
        (stream.tail_coordinates().reference() as isize + offset_shift.reference)
            .try_into()
            .unwrap(),
        query.sequence.chars(),
        stream.stream_iter(),
        true,
        false,
    );

    renderer
        .render_without_names(
            stdout(),
            reference_c
                .as_ref()
                .map(|reference_c| reference_c.label)
                .into_iter()
                .chain([reference.label, query.label])
                .chain(query_c.as_ref().map(|query_c| query_c.label)),
        )
        .unwrap();

    if let Some(reference_c) = reference_c.as_ref() {
        renderer.extend_sequence_with_alignment_and_default_data(
            reference.label,
            reference_c.label,
            (stream.tail_coordinates().reference() as isize + offset_shift.reference)
                .try_into()
                .unwrap(),
            reference_c.sequence.chars(),
            [(
                reference_c.sequence.chars().count(),
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );
    }

    if let Some(query_c) = query_c {
        renderer
            .render_without_names(
                stdout(),
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label, query_c.label]),
            )
            .unwrap();

        renderer.extend_sequence_with_alignment_and_default_data(
            query.label,
            query_c.label,
            (stream.tail_coordinates().query() as isize + offset_shift.query)
                .try_into()
                .unwrap(),
            query_c.sequence.chars(),
            [(
                query_c.sequence.chars().count(),
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );
    }
}

fn render_ts_base<SequenceName: Eq + Ord + Clone>(
    renderer: &mut MultipairAlignmentRenderer<SequenceName, CharacterData>,
    offset_shift: &mut OffsetShift,
    reference: &LabelledSequence<&SequenceName>,
    query: &LabelledSequence<&SequenceName>,
    reference_c: Option<&LabelledSequence<&SequenceName>>,
    query_c: Option<&LabelledSequence<&SequenceName>>,
    stream: &AlignmentStream,
) {
    debug!(
        "Rendering ts base {} from RQ {}/{} to RQ {}/{}",
        stream.stream_iter().next().unwrap().1,
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

    renderer
        .render_without_names(
            stdout(),
            reference_c
                .as_ref()
                .map(|reference_c| reference_c.label)
                .into_iter()
                .chain([reference.label, query.label])
                .chain(query_c.as_ref().map(|query_c| query_c.label)),
        )
        .unwrap();

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
        [(
            secondary.sequence.chars().count(),
            AlignmentType::PrimaryMatch,
        )],
        false,
        false,
    );

    if let Some(secondary_c) = secondary_c {
        renderer.extend_sequence_with_alignment_and_default_data(
            secondary.label,
            secondary_c.label,
            secondary_offset,
            secondary_c.sequence.chars(),
            [(
                secondary_c.sequence.chars().count(),
                AlignmentType::PrimaryMatch,
            )],
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
        trace!("extension_existing_length: {extension_existing_length}");

        renderer
            .render_without_names(
                stdout(),
                reference_c
                    .as_ref()
                    .map(|reference_c| reference_c.label)
                    .into_iter()
                    .chain([reference.label, query.label])
                    .chain(query_c.as_ref().map(|query_c| query_c.label)),
            )
            .unwrap();

        renderer.extend_sequence_with_alignment(
            anti_secondary.label,
            anti_secondary_c.label,
            anti_secondary_offset,
            extension.clone().take(extension_existing_length),
            Default::default,
            Default::default,
            [(
                extension_existing_length.min(extension.clone().count()),
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );
        renderer.extend_sequence(
            anti_secondary_c.label,
            extension.skip(extension_existing_length),
            Default::default,
        );
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
}

fn make_png(output: impl AsRef<std::path::Path>) {
    let svg_in = output.as_ref();
    let png_out = svg_in.with_extension("png");
    info!("Converting SVG to PNG at {png_out:?}");

    let svg = std::fs::read(svg_in).unwrap();
    let svg = usvg::Tree::from_data(&svg, &Default::default()).unwrap();

    let zoom = 20.0;
    let raster_image_size = svg.size();
    let raster_image_width = (raster_image_size.width().ceil() * zoom) as u32;
    let raster_image_height = (raster_image_size.height().ceil() * zoom) as u32;
    info!("PNG size: {raster_image_width}x{raster_image_height}",);

    let mut raster_image = tiny_skia::Pixmap::new(raster_image_width, raster_image_height).unwrap();
    resvg::render(
        &svg,
        Transform::from_scale(zoom, zoom),
        &mut raster_image.as_mut(),
    );
    raster_image.save_png(png_out).unwrap();
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
