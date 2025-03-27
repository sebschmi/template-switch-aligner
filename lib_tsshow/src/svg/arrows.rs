use std::{collections::BTreeMap, fmt::Display};

use log::debug;
use svg::{
    Document,
    node::element::{Definitions, Group, Marker, Path},
};

use crate::{
    plain_text::mutlipair_alignment_renderer::MultipairAlignmentRenderer, svg::font::typewriter,
};

#[derive(Debug, PartialEq)]
pub struct Arrow {
    pub from: ArrowEndpoint,
    pub to: ArrowEndpoint,
    pub style: ArrowStyle,
}

#[derive(Debug, PartialEq)]
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

    /// Extra offset affecting the `y`-coordinate of this endpoint.
    ///
    /// The offset is applied according to [`direction`](Self::direction).
    pub column_offset: f32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ArrowEndpointDirection {
    Forward,
    Backward,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
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
                column_offset: 0.0,
            },
            to: ArrowEndpoint {
                column: to_column,
                row,
                direction: if from_column <= to_column {
                    ArrowEndpointDirection::Backward
                } else {
                    ArrowEndpointDirection::Forward
                },
                column_offset: 0.0,
            },
            style: ArrowStyle::Direct,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_curved(
        from_column: usize,
        from_column_offset: f32,
        from_row: String,
        from_direction: ArrowEndpointDirection,
        to_column: usize,
        to_column_offset: f32,
        to_row: String,
        to_direction: ArrowEndpointDirection,
    ) -> Self {
        Self {
            from: ArrowEndpoint::new(from_column, from_row, from_direction, from_column_offset),
            to: ArrowEndpoint::new(to_column, to_row, to_direction, to_column_offset),
            style: ArrowStyle::Curved,
        }
    }

    pub fn render(&self, rows: &BTreeMap<String, f32>) -> Group {
        match self.style {
            ArrowStyle::Direct => {
                debug!("Drawing direct arrow from {} to {}", self.from, self.to);

                let stroke_width = 0.15;
                let from_x = self.from.column as f32 * typewriter::FONT.character_width
                    + self.from.column_offset * self.from.direction.sign();
                let from_y =
                    *rows.get(&self.from.row).unwrap() - typewriter::FONT.character_height * 0.3;
                let to_x = self.to.column as f32 * typewriter::FONT.character_width
                    + self.to.column_offset * self.to.direction.sign();
                let to_y =
                    *rows.get(&self.to.row).unwrap() - typewriter::FONT.character_height * 0.3;

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
                let from_x = self.from.column as f32 * typewriter::FONT.character_width
                    + self.from.column_offset * self.from.direction.sign();
                let from_y =
                    *rows.get(&self.from.row).unwrap() - typewriter::FONT.character_height * 0.3;
                let to_x = self.to.column as f32 * typewriter::FONT.character_width
                    + self.to.column_offset * self.to.direction.sign();
                let to_y =
                    *rows.get(&self.to.row).unwrap() - typewriter::FONT.character_height * 0.3;

                let (from_x_control, to_x_control) =
                    match (&self.from.direction, &self.to.direction) {
                        (ArrowEndpointDirection::Forward, ArrowEndpointDirection::Forward) => {
                            let x_control =
                                from_x.max(to_x) + 2.0 * typewriter::FONT.character_width;
                            (x_control, x_control)
                        }
                        (ArrowEndpointDirection::Backward, ArrowEndpointDirection::Backward) => {
                            let x_control =
                                from_x.min(to_x) - 2.0 * typewriter::FONT.character_width;
                            (x_control, x_control)
                        }
                        (ArrowEndpointDirection::Forward, ArrowEndpointDirection::Backward)
                        | (ArrowEndpointDirection::Backward, ArrowEndpointDirection::Forward) => {
                            let x_control_delta = ((from_x - to_x).abs() * 0.1)
                                .max(2.0 * typewriter::FONT.character_width);
                            (
                                from_x + self.from.direction.sign() * x_control_delta,
                                to_x + self.to.direction.sign() * x_control_delta,
                            )
                        }
                    };

                Group::new().add(
                    Path::new()
                        .set("marker-end", "url(#arrow_head_red)")
                        .set("stroke-width", stroke_width)
                        .set("stroke", "#CE2029")
                        .set("fill", "none")
                        .set("stroke-dasharray", stroke_width)
                        .set("d", format!("M {from_x},{from_y} C {from_x_control},{from_y} {to_x_control},{to_y} {to_x},{to_y}")),
                )
            }
        }
    }

    pub fn translate_char_gap_columns_to_real_columns<CharacterData>(
        mut self,
        renderer: &MultipairAlignmentRenderer<String, CharacterData>,
    ) -> Self {
        self.from
            .translate_char_gap_columns_to_real_columns(renderer);
        self.to.translate_char_gap_columns_to_real_columns(renderer);
        self
    }
}

impl ArrowEndpoint {
    pub fn new(
        column: usize,
        row: String,
        direction: ArrowEndpointDirection,
        column_offset: f32,
    ) -> Self {
        Self {
            column,
            row,
            direction,
            column_offset,
        }
    }

    pub fn translate_char_gap_columns_to_real_columns<CharacterData>(
        &mut self,
        renderer: &MultipairAlignmentRenderer<String, CharacterData>,
    ) {
        let sequence = renderer.sequence(&self.row);
        self.column = match self.direction {
            ArrowEndpointDirection::Forward => {
                if self.column == 0 {
                    0
                } else {
                    sequence
                        .translate_offset_without_blanks(self.column - 1)
                        .unwrap()
                        + 1
                }
            }

            ArrowEndpointDirection::Backward => sequence
                .translate_offset_without_blanks(self.column)
                .unwrap(),
        };
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
        Definitions::new()
            .add(
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
            )
            .add(
                Marker::new()
                    .set("id", "arrow_head_red")
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
                            .set("stroke", "#CE2029"),
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
