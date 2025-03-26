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
    Curved,
}

impl Arrow {
    pub fn new_skip(from_column: usize, to_column: usize, row: String) -> Self {
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

    pub fn new_curved(
        from_column: usize,
        from_row: String,
        from_direction: ArrowEndpointDirection,
        to_column: usize,
        to_row: String,
        to_direction: ArrowEndpointDirection,
    ) -> Self {
        Self {
            from: ArrowEndpoint::new(from_column, from_row, from_direction),
            to: ArrowEndpoint::new(to_column, to_row, to_direction),
            style: ArrowStyle::Curved,
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
                        .set("stroke-dasharray", stroke_width)
                        .set("d", format!("M {from_x},{from_y} L {to_x},{to_y}")),
                )
            }

            ArrowStyle::Curved => {
                debug!("Drawing curved arrow from {} to {}", self.from, self.to);

                let stroke_width = 0.15;
                let from_x = self.from.column as f32 * CHARACTER_WIDTH;
                let from_y = *rows.get(&self.from.row).unwrap() - CHARACTER_HEIGHT * 0.3;
                let to_x = self.to.column as f32 * CHARACTER_WIDTH - 2.0 * stroke_width;
                let to_y = *rows.get(&self.to.row).unwrap() - CHARACTER_HEIGHT * 0.3;

                let from_x_control = from_x + 5.0 * self.from.direction.sign();
                let to_x_control = to_x + 5.0 * self.to.direction.sign();

                Group::new().add(
                    Path::new()
                        .set("marker-end", "url(#arrow_head)")
                        .set("stroke-width", stroke_width)
                        .set("stroke", "black")
                        .set("fill", "none")
                        .set("stroke-dasharray", stroke_width)
                        .set("d", format!("M {from_x},{from_y} C {from_x_control},{from_y} {to_x_control},{to_y} {to_x},{to_y}")),
                )
            }
        }
    }
}

impl ArrowEndpoint {
    pub fn new(column: usize, row: String, direction: ArrowEndpointDirection) -> Self {
        Self {
            column,
            row,
            direction,
        }
    }
}

impl ArrowEndpointDirection {
    /// Returns `1.0` if the direction is forward, and `-1.0` if the direction is backward.
    pub fn sign(&self) -> f32 {
        match self {
            ArrowEndpointDirection::Forward => 1.0,
            ArrowEndpointDirection::Backward => -1.0,
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

impl Display for Arrow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}({} -> {})",
            match self.style {
                ArrowStyle::Direct => "D",
                ArrowStyle::Curved => "C",
            },
            self.from,
            self.to
        )
    }
}

impl Display for ArrowEndpoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}/{}", self.column, self.row, self.direction)
    }
}

impl Display for ArrowEndpointDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ArrowEndpointDirection::Forward => "+",
                ArrowEndpointDirection::Backward => "-",
            }
        )
    }
}
