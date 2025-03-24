use std::{collections::BTreeMap, fmt::Display};

use log::debug;
use svg::{
    Document,
    node::element::{Definitions, Group, Marker, Path},
};

use super::font::{CHARACTER_HEIGHT, CHARACTER_WIDTH};

pub struct Arrow {
    pub from: ArrowEndpoint,
    pub to: ArrowEndpoint,
    pub style: ArrowStyle,
}

pub struct ArrowEndpoint {
    /// The column the endpoint is at, which will be translated to an `x`-coordinate.
    ///
    /// The column is interpreted as the left side of the bounding box of a character.
    pub column: usize,

    /// The row the endpoint is at, which will be translated to a `y`-coordinate.
    ///
    /// The row is interpreted as the horizontal center of the bounding box of a character.
    pub row: String,

    /// The direction from which the arrow travels to this endpoint.
    ///
    /// **Note:** this has nothing to do with where the head and tail of the arrow are,
    /// but it decides where the curve the arrow is made of goes.
    pub direction: ArrowEndpointDirection,
}

pub enum ArrowEndpointDirection {
    Forward,
    Backward,
}

pub enum ArrowStyle {
    Direct,
}

impl Arrow {
    pub fn new_skip(from_column: usize, to_column: usize, row: String) -> Arrow {
        Self {
            from: ArrowEndpoint {
                column: from_column,
                row: row.clone(),
                direction: if from_column <= to_column {
                    ArrowEndpointDirection::Forward
                } else {
                    ArrowEndpointDirection::Backward
                },
            },
            to: ArrowEndpoint {
                column: to_column,
                row,
                direction: if from_column <= to_column {
                    ArrowEndpointDirection::Backward
                } else {
                    ArrowEndpointDirection::Forward
                },
            },
            style: ArrowStyle::Direct,
        }
    }

    pub fn render(&self, rows: &BTreeMap<String, f32>) -> Group {
        match self.style {
            ArrowStyle::Direct => {
                debug!("Drawing direct arrow from {} to {}", self.from, self.to);

                let stroke_width = 0.15;
                let from_x = self.from.column as f32 * CHARACTER_WIDTH;
                let from_y = *rows.get(&self.from.row).unwrap() - CHARACTER_HEIGHT * 0.3;
                let to_x = self.to.column as f32 * CHARACTER_WIDTH - 2.0 * stroke_width;
                let to_y = *rows.get(&self.to.row).unwrap() - CHARACTER_HEIGHT * 0.3;

                Group::new().add(
                    Path::new()
                        .set("marker-end", "url(#arrow_head)")
                        .set("stroke-width", stroke_width)
                        .set("stroke", "black")
                        .set("fill", "none")
                        .set("stroke-dasharray", format!("{}", 1.0 / 1.5))
                        .set("d", format!("M {from_x},{from_y} L {to_x},{to_y}")),
                )
            }
        }
    }
}

pub fn add_arrow_defs(svg: Document) -> Document {
    svg.add(
        Definitions::new().add(
            Marker::new()
                .set("id", "arrow_head")
                .set("viewBox", "0 0 10 10")
                .set("orient", "auto-start-reverse")
                .set("markerWidth", 10)
                .set("markerHeight", 10)
                .set("refX", 10)
                .set("refY", 5)
                .add(
                    Path::new()
                        .set("d", "M 1 1 L 10 5 L 1 9")
                        .set("fill", "none")
                        .set("stroke", "black"),
                ),
        ),
    )
}

impl Display for ArrowEndpoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}", self.column, self.row)
    }
}
