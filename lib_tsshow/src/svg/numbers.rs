use std::fmt::Display;

use svg::node::element::{Group, Path};

use crate::{
    plain_text::mutlipair_alignment_renderer::{Character, NoCharacterData},
    svg::font::sans_serif_mono,
    ts_arrangement::index_types::ArrangementColumn,
};

use super::{
    SvgLocation,
    font::{svg_string, typewriter},
};

const BORDER_STROKE_WIDTH: f32 = 0.2;
const BORDER_PADDING: f32 = 0.4;
const BORDER_RADIUS: f32 = 1.2;

pub struct Number<Row> {
    number: String,
    column: ArrangementColumn,
    row: Row,
    align: NumberAlignment,
    scale: f32,
}

pub enum NumberAlignment {
    Left,
    Right,
}

impl<Row> Number<Row> {
    pub fn new(
        number: impl ToString,
        column: ArrangementColumn,
        row: Row,
        align: NumberAlignment,
        scale: f32,
    ) -> Self {
        Self {
            number: number.to_string(),
            column,
            row,
            align,
            scale,
        }
    }

    pub fn render(&self, mut rows: impl FnMut(&Row) -> f32) -> Group {
        let x = self.column.primitive() as f32 * typewriter::FONT.character_width
            + (BORDER_PADDING + 0.5 * BORDER_STROKE_WIDTH) * self.scale
            - match self.align {
                NumberAlignment::Left => 0.0,
                NumberAlignment::Right => self.width(),
            };
        let y = rows(&self.row);
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        let label_height = sans_serif_mono::FONT.character_height;

        let border_x = BORDER_RADIUS - BORDER_PADDING;
        let border_y = -BORDER_PADDING - label_height;
        let border_h = label_width + 2.0 * (BORDER_PADDING - BORDER_RADIUS);
        let border_v = label_height + 2.0 * (BORDER_PADDING - BORDER_RADIUS);

        let border_d = format!(
            "M {border_x},{border_y} a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 -{BORDER_RADIUS},{BORDER_RADIUS} v {border_v}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 {BORDER_RADIUS},{BORDER_RADIUS} h {border_h}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 {BORDER_RADIUS},-{BORDER_RADIUS} v -{border_v}\
             a {BORDER_RADIUS},{BORDER_RADIUS} 90 0 0 -{BORDER_RADIUS},-{BORDER_RADIUS} h -{border_h}",
        );

        Group::new()
            .set(
                "transform",
                format!("translate({x},{y}) scale({})", self.scale),
            )
            .add(
                Path::new()
                    .set("fill", "none")
                    .set("stroke", "black")
                    .set("stroke-width", BORDER_STROKE_WIDTH)
                    .set("d", border_d),
            )
            .add(
                Group::new().add(svg_string(
                    self.number
                        .chars()
                        .map(|c| Character::new_char(c, NoCharacterData)),
                    &SvgLocation {
                        x: 0.0,
                        y: sans_serif_mono::FONT.character_height * -0.12,
                    },
                    &sans_serif_mono::FONT,
                )),
            )
    }

    pub fn width(&self) -> f32 {
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        (label_width + BORDER_STROKE_WIDTH + 2.0 * BORDER_PADDING) * self.scale
    }
}

impl<Row: Display> Display for Number<Row> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} at {}/{}/{}",
            self.number, self.column, self.row, self.align
        )
    }
}

impl Display for NumberAlignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                NumberAlignment::Left => "L",
                NumberAlignment::Right => "R",
            }
        )
    }
}
