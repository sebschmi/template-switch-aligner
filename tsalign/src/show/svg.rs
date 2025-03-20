use std::iter;

use font::{CHARACTER_HEIGHT, CHARACTER_WIDTH, svg_string};
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::AlignmentResult,
        template_switch_distance::{AlignmentType, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info, trace, warn};
use resvg::{
    tiny_skia,
    usvg::{self, Transform},
};
use svg::{Document, node::element::Group};

use crate::show::{
    alignment_stream::AlignmentStream, mutlipair_alignment_renderer::MultipairAlignmentRenderer,
};

mod font;

struct SvgLocation {
    pub x: f32,
    pub y: f32,
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

    if !alignment.iter().any(|(_, alignment_type)| {
        matches!(alignment_type, AlignmentType::TemplateSwitchEntrance { .. })
    }) {
        warn!("No template switches found");
    }

    let reference_label = "Reference".to_string();
    let query_label = "Query".to_string();

    let reference = &result.statistics().sequences.reference;
    let query = &result.statistics().sequences.query;

    debug!("Rendering reference and query");
    let mut stream = AlignmentStream::new();
    let mut offset_stream = AlignmentStream::new();
    let mut renderer = MultipairAlignmentRenderer::new(reference_label.clone(), "");
    renderer.add_independent_sequence(query_label.clone(), "");
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
                &reference_label,
                &query_label,
                reference,
                query,
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

            render_ts_base(
                &mut renderer,
                &mut reference_render_offset_shift,
                &mut query_render_offset_shift,
                &reference_label,
                &query_label,
                reference,
                query,
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
            &reference_label,
            &query_label,
            reference,
            query,
            &stream,
        );
    }

    let reference = svg_string(
        renderer.sequence(&reference_label).characters(),
        &SvgLocation { x: 0.0, y: 0.0 },
    );
    let query = svg_string(
        renderer.sequence(&query_label).characters(),
        &SvgLocation {
            x: 0.0,
            y: CHARACTER_HEIGHT,
        },
    );
    let rq_group = Group::new()
        .set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform())
        .add(reference)
        .add(query);

    let mut view_box_width = offset_stream.len() as f32 * CHARACTER_WIDTH + 20.0;
    let mut view_box_height = 3.0 * CHARACTER_HEIGHT + 20.0;

    let mut svg = Document::new().add(rq_group);

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
        let mut renderer = MultipairAlignmentRenderer::new(
            reference_label.clone(),
            &no_ts_result.statistics().sequences.reference,
        );

        debug!("Adding query");
        renderer.add_aligned_sequence(
            &reference_label,
            0,
            query_label.clone(),
            &no_ts_result.statistics().sequences.query,
            no_ts_alignment.iter().copied(),
            true,
            false,
        );

        debug!("Rendering");
        let reference = svg_string(
            renderer.sequence(&reference_label).characters(),
            &SvgLocation { x: 0.0, y: 0.0 },
        );
        let query = svg_string(
            renderer.sequence(&query_label).characters(),
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
    renderer: &mut MultipairAlignmentRenderer<SequenceName>,
    reference_render_offset_shift: isize,
    reference_label: &SequenceName,
    query_label: &SequenceName,
    reference: &str,
    query: &str,
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

    let reference =
        &reference[stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
    let query = &query[stream.tail_coordinates().query()..stream.head_coordinates().query()];

    renderer.extend_sequence(reference_label, reference);
    renderer.extend_sequence_with_alignment(
        reference_label,
        query_label,
        (stream.tail_coordinates().reference() as isize + reference_render_offset_shift)
            .try_into()
            .unwrap(),
        query.chars(),
        stream.stream_iter(),
        true,
        false,
    );
}

#[expect(clippy::too_many_arguments)]
fn render_ts_base<SequenceName: Eq + Ord>(
    renderer: &mut MultipairAlignmentRenderer<SequenceName>,
    reference_render_offset_shift: &mut isize,
    query_render_offset_shift: &mut isize,
    reference_label: &SequenceName,
    query_label: &SequenceName,
    reference: &str,
    query: &str,
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

    let reference =
        &reference[stream.tail_coordinates().reference()..stream.head_coordinates().reference()];
    let query = &query[stream.tail_coordinates().query()..stream.head_coordinates().query()];

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
            renderer.extend_sequence(reference_label, reference);
            renderer.extend_sequence_with_alignment(
                reference_label,
                query_label,
                (stream.tail_coordinates().reference() as isize + *reference_render_offset_shift)
                    .try_into()
                    .unwrap(),
                iter::repeat_n(' ', reference.chars().count()),
                [(reference.chars().count(), AlignmentType::PrimaryMatch)],
                false,
                false,
            );
            *query_render_offset_shift += length_difference;
        }
        TemplateSwitchSecondary::Query => {
            renderer.extend_sequence(query_label, query);
            renderer.extend_sequence_with_alignment(
                query_label,
                reference_label,
                (stream.tail_coordinates().query() as isize + *query_render_offset_shift)
                    .try_into()
                    .unwrap(),
                iter::repeat_n(' ', query.chars().count()),
                [(query.chars().count(), AlignmentType::PrimaryMatch)],
                false,
                false,
            );
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
