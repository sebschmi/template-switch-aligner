use std::iter;

use font::{CHARACTER_HEIGHT, CHARACTER_WIDTH, CharacterData, svg_string};
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
    Document,
    node::element::{Circle, Group},
};

use crate::show::{
    alignment_stream::AlignmentStream,
    mutlipair_alignment_renderer::{Character, MultipairAlignmentRenderer},
};

use super::alignment_stream::AlignmentCoordinates;

mod font;

struct SvgLocation {
    pub x: f32,
    pub y: f32,
}

struct LabelledSequence<'a, Label> {
    label: Label,
    sequence: &'a str,
}

struct TemplateSwitch {
    tail: AlignmentCoordinates,
    head: AlignmentCoordinates,
    alignment: Vec<(usize, AlignmentType)>,
}

impl SvgLocation {
    pub fn as_transform(&self) -> String {
        format!("translate({} {})", self.x, self.y)
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
            }
            (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Query) => {
                has_secondary_query_ts = true
            }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Reference) => {
                has_secondary_reference_ts = true
            }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Query) => {
                has_self_query_ts = true
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

    let mut renderer = MultipairAlignmentRenderer::new_empty();
    renderer.add_empty_independent_sequence(reference_label.clone());
    renderer.add_empty_independent_sequence(query_label.clone());
    if has_secondary_reference_ts {
        renderer.add_empty_independent_sequence(reference_c_label.clone());
    }
    if has_secondary_query_ts {
        renderer.add_empty_independent_sequence(query_c_label.clone());
    }

    let mut reference_render_offset_shift = 0;
    let mut query_render_offset_shift = 0;

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
        trace!("Outer {alignment_type}");

        if matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. }) {
            render_inter_ts(
                &mut renderer,
                reference_render_offset_shift,
                query_render_offset_shift,
                &reference,
                &query,
                has_secondary_reference_ts.then_some(&reference_c),
                has_secondary_query_ts.then_some(&query_c),
                &stream,
            );

            stream.clear();

            // Accumulate TS
            for alignment_type in alignment_iter.by_ref() {
                trace!("Skipping {alignment_type}");

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
                &mut reference_render_offset_shift,
                &mut query_render_offset_shift,
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
            reference_render_offset_shift,
            query_render_offset_shift,
            &reference,
            &query,
            has_secondary_reference_ts.then_some(&reference_c),
            has_secondary_query_ts.then_some(&query_c),
            &stream,
        );
    }

    debug!("Rendering TSes");
    let mut ts_rq_labels = Vec::new();
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

        match (primary, secondary) {
            (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Reference) => todo!(),
            (TemplateSwitchPrimary::Reference, TemplateSwitchSecondary::Query) => {
                let label = format!("{reference_label} TS{}", index + 1);
                let inner_sequence = &reference.sequence[tail.reference()..head.reference()];
                let inner_sequence_r: String = inner_sequence.chars().rev().collect();
                let inner_offset: usize = (tail.reference() as isize + first_offset)
                    .try_into()
                    .unwrap();

                renderer.add_aligned_sequence(
                    query_c.label,
                    inner_offset,
                    label.clone(),
                    inner_sequence_r
                        .chars()
                        .map(|c| Character::new_char(c, CharacterData::default())),
                    CharacterData::default,
                    CharacterData::default,
                    alignment[1..alignment.len() - 1].iter().copied(),
                    true,
                    false,
                );

                ts_rq_labels.push(label);
            }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Reference) => { /* TODO */ }
            (TemplateSwitchPrimary::Query, TemplateSwitchSecondary::Query) => todo!(),
        }
    }

    debug!("Rendering SVG");
    let mut rq_group =
        Group::new().set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform());
    let mut y = 0.0;

    if has_secondary_reference_ts {
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

        for ts_label in &ts_rq_labels {
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

#[expect(clippy::too_many_arguments)]
fn render_inter_ts<SequenceName: Eq + Ord>(
    renderer: &mut MultipairAlignmentRenderer<SequenceName, CharacterData>,
    reference_render_offset_shift: isize,
    query_render_offset_shift: isize,
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
    debug!("Reference render offset shift: {reference_render_offset_shift}");

    let reference_sequence = &reference.sequence
        [stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
    let query_sequence =
        &query.sequence[stream.tail_coordinates().query()..stream.head_coordinates().query()];

    renderer.extend_sequence_with_default_data(reference.label, reference_sequence.chars());
    renderer.extend_sequence_with_alignment_and_default_data(
        reference.label,
        query.label,
        (stream.tail_coordinates().reference() as isize + reference_render_offset_shift)
            .try_into()
            .unwrap(),
        query_sequence.chars(),
        stream.stream_iter(),
        true,
        false,
    );

    if let Some(reference_c) = reference_c {
        let reference_c_sequence = &reference_c.sequence
            [stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
        renderer.extend_sequence_with_alignment_and_default_data(
            reference.label,
            reference_c.label,
            (stream.tail_coordinates().reference() as isize + reference_render_offset_shift)
                .try_into()
                .unwrap(),
            reference_c_sequence.chars(),
            [(
                reference_c_sequence.chars().count(),
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );
    }

    if let Some(query_c) = query_c {
        let query_c_sequence =
            &query_c.sequence[stream.tail_coordinates().query()..stream.head_coordinates().query()];
        renderer.extend_sequence_with_alignment_and_default_data(
            query.label,
            query_c.label,
            (stream.tail_coordinates().query() as isize + query_render_offset_shift)
                .try_into()
                .unwrap(),
            query_c_sequence.chars(),
            [(
                query_c_sequence.chars().count(),
                AlignmentType::PrimaryMatch,
            )],
            false,
            false,
        );
    }
}

#[expect(clippy::too_many_arguments)]
fn render_ts_base<SequenceName: Eq + Ord>(
    renderer: &mut MultipairAlignmentRenderer<SequenceName, CharacterData>,
    reference_render_offset_shift: &mut isize,
    query_render_offset_shift: &mut isize,
    reference: &LabelledSequence<&SequenceName>,
    query: &LabelledSequence<&SequenceName>,
    reference_c: Option<&LabelledSequence<&SequenceName>>,
    query_c: Option<&LabelledSequence<&SequenceName>>,
    stream: &AlignmentStream,
) {
    debug!(
        "Rendering ts base from RQ {}/{} to RQ {}/{}",
        stream.tail_coordinates().reference(),
        stream.tail_coordinates().query(),
        stream.head_coordinates().reference(),
        stream.head_coordinates().query()
    );
    debug!("Render offset shifts RQ {reference_render_offset_shift}/{query_render_offset_shift}");

    let reference_sequence = &reference.sequence
        [stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
    let query_sequence =
        &query.sequence[stream.tail_coordinates().query()..stream.head_coordinates().query()];

    let mut alignment = stream.stream_iter();
    let (_, AlignmentType::TemplateSwitchEntrance { secondary, .. }) = alignment.next().unwrap()
    else {
        panic!("Wrong template switch entrance");
    };
    let (_, AlignmentType::TemplateSwitchExit { length_difference }) = alignment.last().unwrap()
    else {
        panic!("Wrong template switch exit")
    };

    match secondary {
        TemplateSwitchSecondary::Reference => {
            renderer.extend_sequence_with_default_data(reference.label, reference_sequence.chars());
            renderer.extend_sequence_with_alignment_and_default_data(
                reference.label,
                query.label,
                (stream.tail_coordinates().reference() as isize + *reference_render_offset_shift)
                    .try_into()
                    .unwrap(),
                iter::repeat_n(' ', reference_sequence.chars().count()),
                [(
                    reference_sequence.chars().count(),
                    AlignmentType::PrimaryMatch,
                )],
                false,
                false,
            );

            if let Some(reference_c) = reference_c {
                let reference_c_sequence = &reference_c.sequence
                    [stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
                renderer.extend_sequence_with_alignment_and_default_data(
                    reference.label,
                    reference_c.label,
                    (stream.tail_coordinates().reference() as isize
                        + *reference_render_offset_shift)
                        .try_into()
                        .unwrap(),
                    reference_c_sequence.chars(),
                    [(
                        reference_c_sequence.chars().count(),
                        AlignmentType::PrimaryMatch,
                    )],
                    false,
                    false,
                );
            }
            if let Some(query_c) = query_c {
                renderer.extend_sequence_with_alignment_and_default_data(
                    query.label,
                    query_c.label,
                    (stream.tail_coordinates().query() as isize + *query_render_offset_shift)
                        .try_into()
                        .unwrap(),
                    iter::repeat_n(' ', reference_sequence.chars().count()),
                    [(
                        reference_sequence.chars().count(),
                        AlignmentType::PrimaryMatch,
                    )],
                    false,
                    false,
                );
            }

            *query_render_offset_shift += length_difference;
        }
        TemplateSwitchSecondary::Query => {
            renderer.extend_sequence_with_default_data(query.label, query_sequence.chars());
            renderer.extend_sequence_with_alignment_and_default_data(
                query.label,
                reference.label,
                (stream.tail_coordinates().query() as isize + *query_render_offset_shift)
                    .try_into()
                    .unwrap(),
                iter::repeat_n(' ', query_sequence.chars().count()),
                [(query_sequence.chars().count(), AlignmentType::PrimaryMatch)],
                false,
                false,
            );

            if let Some(reference_c) = reference_c {
                renderer.extend_sequence_with_alignment_and_default_data(
                    reference.label,
                    reference_c.label,
                    (stream.tail_coordinates().reference() as isize
                        + *reference_render_offset_shift)
                        .try_into()
                        .unwrap(),
                    iter::repeat_n(' ', query_sequence.chars().count()),
                    [(query_sequence.chars().count(), AlignmentType::PrimaryMatch)],
                    false,
                    false,
                );
            }
            if let Some(query_c) = query_c {
                let query_c_sequence = &query_c.sequence
                    [stream.tail_coordinates().query()..stream.head_coordinates().query()];
                renderer.extend_sequence_with_alignment_and_default_data(
                    query.label,
                    query_c.label,
                    (stream.tail_coordinates().query() as isize + *query_render_offset_shift)
                        .try_into()
                        .unwrap(),
                    query_c_sequence.chars(),
                    [(
                        query_c_sequence.chars().count(),
                        AlignmentType::PrimaryMatch,
                    )],
                    false,
                    false,
                );
            }

            *reference_render_offset_shift += length_difference;
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
