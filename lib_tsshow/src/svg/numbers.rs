use std::collections::BTreeMap;

use svg::node::element::{Group, Path};

use crate::{
    plain_text::mutlipair_alignment_renderer::{Character, NoCharacterData},
    svg::font::sans_serif_mono,
};

use super::{SvgLocation, font::svg_string};

const BORDER_STROKE_WIDTH: f32 = 0.1;
const BORDER_PADDING: f32 = 0.2;

pub struct Number {
    number: String,
    column: usize,
    row: String,
    align: NumberAlign,
}

pub enum NumberAlign {
    Left,
    Right,
}

impl Number {
    pub fn new(number: impl ToString, column: usize, row: String, align: NumberAlign) -> Self {
        Self {
            number: number.to_string(),
            column,
            row,
            align,
        }
    }

    pub fn render(&self, rows: &BTreeMap<String, f32>) -> Group {
        let x = self.column as f32 * sans_serif_mono::FONT.character_width;
        let y = *rows.get(&self.row).unwrap();
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        let label_height = sans_serif_mono::FONT.character_height;

        Group::new()
            .set("transform", SvgLocation {x, y}.as_transform())
            .add(
            Path::new()
                    .set("fill", "none")
                    .set("stroke", "black")
                    .set("stroke-width", BORDER_STROKE_WIDTH)
                    .set("d", format!("M {BORDER_PADDING},0 a {BORDER_PADDING},{BORDER_PADDING} 90 0 0 -{BORDER_PADDING},{BORDER_PADDING} v {label_height}\
                        a {BORDER_PADDING},{BORDER_PADDING} 90 0 0 {BORDER_PADDING},{BORDER_PADDING} h {label_width}\
                        a {BORDER_PADDING},{BORDER_PADDING} 90 0 0 {BORDER_PADDING},-{BORDER_PADDING} v -{label_height}\
                        a {BORDER_PADDING},{BORDER_PADDING} 90 0 0 -{BORDER_PADDING},-{BORDER_PADDING} h -{label_width}"))
            .add(svg_string(
                self.number.chars().map(|c| Character::new_char(c, NoCharacterData)),
                &SvgLocation { x: BORDER_PADDING, y: BORDER_PADDING }, 
                &sans_serif_mono::FONT
            )),
        )
    }

    pub fn width(&self) -> f32 {
        let label_width =
            self.number.chars().count() as f32 * sans_serif_mono::FONT.character_width;
        label_width + BORDER_STROKE_WIDTH + 2.0 * BORDER_PADDING
    }
}
