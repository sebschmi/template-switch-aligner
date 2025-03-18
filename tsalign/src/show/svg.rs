use font::{CHARACTER_HEIGHT, CHARACTER_WIDTH, svg_string};
use lib_tsalign::{
    a_star_aligner::{alignment_result::AlignmentResult, template_switch_distance::AlignmentType},
    costs::U64Cost,
};
use log::{debug, info, warn};
use svg::Document;

use crate::show::mutlipair_alignment_renderer::MultipairAlignmentRenderer;

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
) {
    let output = output.as_ref();
    info!("Creating template switch SVG at {output:?}");

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

        let reference_label = "Reference".to_string();
        let query_label = "Query".to_string();

        debug!("Creating no-ts renderer");
        let mut renderer = MultipairAlignmentRenderer::new(
            reference_label.clone(),
            &no_ts_result.statistics().sequences.reference,
        );

        debug!("Adding primary");
        renderer.add_aligned_sequence(
            &reference_label,
            0,
            query_label.clone(),
            &no_ts_result.statistics().sequences.query,
            no_ts_alignment,
            true,
            false,
        );

        debug!("Rendering");
        let reference = svg_string(
            renderer.sequence(&reference_label).characters(),
            &SvgLocation {
                x: 10.0,
                y: 1.0 * CHARACTER_HEIGHT,
            },
        );
        let query = svg_string(
            renderer.sequence(&query_label).characters(),
            &SvgLocation {
                x: 10.0,
                y: 2.0 * CHARACTER_HEIGHT,
            },
        );

        let svg = Document::new()
            .set(
                "viewBox",
                (
                    0,
                    0,
                    renderer.sequence(&reference_label).len() as f32 * CHARACTER_WIDTH + 20.0,
                    3.0 * CHARACTER_HEIGHT,
                ),
            )
            .add(reference)
            .add(query);

        svg::save(output, &svg).unwrap();
    } else {
        debug!("No no-ts alignment given, skipping");
    }
}
