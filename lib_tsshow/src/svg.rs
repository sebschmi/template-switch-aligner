use std::{collections::HashMap, io::Write};

use arrows::{Arrow, ArrowEndpointDirection, add_arrow_defs};
use font::{CharacterData, svg_string, typewriter};
use indexed_str::IndexedStr;
use lib_tsalign::{
    a_star_aligner::{
        alignment_result::AlignmentResult,
        template_switch_distance::{AlignmentType, TemplateSwitchPrimary, TemplateSwitchSecondary},
    },
    costs::U64Cost,
};
use log::{debug, info};
use numbers::{Number, NumberAlignment};
use svg::{
    Document,
    node::element::{Circle, Group, Text},
};

use crate::{
    error::{Error, Result},
    plain_text::mutlipair_alignment_renderer::Character,
    ts_arrangement::{
        TsArrangement,
        complement::ComplementChar,
        index_types::ArrangementColumn,
        inner::InnerChar,
        row::TsArrangementRow,
        source::{SourceChar, TsSourceArrangement},
        template_switch::TemplateSwitch,
    },
};

mod arrows;
mod font;
mod indexed_str;
pub mod labelled_sequence;
mod numbers;
mod offset_shift;

const COPY_COLORS: &[&str] = &["#00CC00", "#009900", "#006600", "#003300"];
const COMPLEMENT_SOURCE_HIDDEN_COLOR: &str = "grey";
const TS_RUNNING_NUMBER: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

struct SvgLocation {
    pub x: f32,
    pub y: f32,
}

impl SvgLocation {
    pub fn as_transform(&self) -> String {
        format!("translate({} {})", self.x, self.y)
    }
}

pub struct SvgConfig {
    pub render_arrows: bool,
    pub render_more_complement: bool,
}

pub fn create_ts_svg(
    output: impl Write,
    result: &AlignmentResult<AlignmentType, U64Cost>,
    no_ts_result: &Option<AlignmentResult<AlignmentType, U64Cost>>,
    config: &SvgConfig,
) -> Result<()> {
    info!("Creating template switch SVG");

    let AlignmentResult::WithTarget {
        alignment,
        statistics,
    } = result
    else {
        return Err(Error::AlignmentHasNoTarget);
    };
    debug!("Alignment: {alignment:?}");

    let reference = &statistics.sequences.reference;
    let query = &statistics.sequences.query;
    let reference_c: String = statistics.sequences.reference_rc.chars().rev().collect();
    let query_c: String = statistics.sequences.query_rc.chars().rev().collect();

    debug!("Computing TS arrangement");
    let mut ts_arrangement =
        TsArrangement::new(reference.len(), query.len(), alignment.iter_flat_cloned());
    ts_arrangement.remove_empty_columns();

    if config.render_more_complement {
        ts_arrangement.show_complete_complements_if_used();
    }

    let ts_arrangement = ts_arrangement;

    let no_ts_arrangement = no_ts_result
        .as_ref()
        .map(|no_ts_result| {
            debug!("Computing no-TS arrangement");
            let AlignmentResult::WithTarget {
                alignment,
                statistics,
            } = no_ts_result
            else {
                return Err(Error::NoTsAlignmentHasNoTarget);
            };

            let reference = &statistics.sequences.reference;
            let query = &statistics.sequences.query;
            Ok(TsSourceArrangement::new(
                reference.len(),
                query.len(),
                alignment.iter_flat_cloned(),
                &mut Vec::new(),
            ))
        })
        .transpose()?;

    debug!("Computing TS arrows and labels");
    let mut arrows = Vec::new();
    let mut numbers = Vec::new();

    for (
        inner_identifier,
        TemplateSwitch {
            index: ts_index,
            primary,
            secondary,
            sp1_reference,
            sp1_query,
            sp4_reference,
            sp4_query,
            sp2_secondary,
            sp3_secondary,
            ..
        },
    ) in ts_arrangement.template_switches()
    {
        let (anti_primary_sp4_minus_one, anti_primary_sp4, anti_primary_row) = match primary {
            TemplateSwitchPrimary::Reference => (
                sp4_query.checked_sub(1).map(|column| {
                    ts_arrangement.query_arrangement_char_to_arrangement_column(column) + 1usize
                }),
                ts_arrangement.query_arrangement_char_to_arrangement_column(*sp4_query),
                TsArrangementRow::Query,
            ),
            TemplateSwitchPrimary::Query => (
                sp4_reference.checked_sub(1).map(|column| {
                    ts_arrangement.reference_arrangement_char_to_arrangement_column(column) + 1usize
                }),
                ts_arrangement.reference_arrangement_char_to_arrangement_column(*sp4_reference),
                TsArrangementRow::Reference,
            ),
        };

        if let Some(anti_primary_sp4_minus_one) = anti_primary_sp4_minus_one {
            if (anti_primary_sp4 - anti_primary_sp4_minus_one) > ArrangementColumn::ZERO {
                arrows.push(Arrow::new_skip(
                    anti_primary_sp4_minus_one,
                    anti_primary_sp4,
                    anti_primary_row,
                ));
            }
        }

        let (primary_sp1, primary_sp4, primary_row) = match primary {
            TemplateSwitchPrimary::Reference => (
                ts_arrangement.reference_arrangement_char_to_arrangement_column(*sp1_reference),
                ts_arrangement.reference_arrangement_char_to_arrangement_column(*sp4_reference),
                TsArrangementRow::Reference,
            ),
            TemplateSwitchPrimary::Query => (
                ts_arrangement.query_arrangement_char_to_arrangement_column(*sp1_query),
                ts_arrangement.query_arrangement_char_to_arrangement_column(*sp4_query),
                TsArrangementRow::Query,
            ),
        };

        let (secondary_sp2, secondary_sp3, inner_row) = match secondary {
            TemplateSwitchSecondary::Reference => (
                ts_arrangement.reference_source_to_arrangement_column(*sp2_secondary),
                ts_arrangement.reference_source_to_arrangement_column(*sp3_secondary),
                TsArrangementRow::ReferenceInner {
                    index: inner_identifier,
                },
            ),
            TemplateSwitchSecondary::Query => (
                ts_arrangement.query_source_to_arrangement_column(*sp2_secondary),
                ts_arrangement.query_source_to_arrangement_column(*sp3_secondary),
                TsArrangementRow::QueryInner {
                    index: inner_identifier,
                },
            ),
        };

        let running_number = TS_RUNNING_NUMBER.chars().nth(*ts_index).unwrap();
        let number1 = Number::new(
            format!("{running_number}1"),
            primary_sp1,
            primary_row,
            NumberAlignment::Left,
            0.5,
        );
        let number2 = Number::new(
            format!("{running_number}2"),
            secondary_sp2,
            inner_row,
            NumberAlignment::Left,
            0.5,
        );
        let number3 = Number::new(
            format!("{running_number}3"),
            secondary_sp3,
            inner_row,
            NumberAlignment::Right,
            0.5,
        );
        let number4 = Number::new(
            format!("{running_number}4"),
            primary_sp4,
            primary_row,
            NumberAlignment::Right,
            0.5,
        );

        if config.render_arrows {
            arrows.push(Arrow::new_curved(
                primary_sp1,
                number1.width(),
                primary_row,
                ArrowEndpointDirection::Forward,
                secondary_sp2,
                number2.width(),
                inner_row,
                ArrowEndpointDirection::Forward,
            ));

            arrows.push(Arrow::new_curved(
                secondary_sp3,
                number3.width(),
                inner_row,
                ArrowEndpointDirection::Backward,
                primary_sp4,
                number4.width(),
                primary_row,
                ArrowEndpointDirection::Backward,
            ));
        }

        numbers.extend([number1, number2, number3, number4]);
    }

    debug!("Creating SVG root");
    let mut svg = Document::new();
    svg = add_arrow_defs(svg);
    svg = svg.add(Circle::new().set("r", 1e5).set("fill", "white"));

    debug!("Rendering TS arrangement");
    let mut ts_group =
        Group::new().set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform());

    let reference = IndexedStr::new(reference);
    let query = IndexedStr::new(query);
    let reference_c = IndexedStr::new(&reference_c);
    let query_c = IndexedStr::new(&query_c);

    let mut rows = HashMap::new();
    let mut y = 0.0;
    for (identifier, reference_inner) in ts_arrangement.reference_inners().iter().rev() {
        ts_group = ts_group.add(svg_string(
            reference_inner
                .sequence()
                .iter_values()
                .map(render_inner_char(
                    match reference_inner.template_switch().primary {
                        TemplateSwitchPrimary::Reference => &reference,
                        TemplateSwitchPrimary::Query => &query,
                    },
                )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(TsArrangementRow::ReferenceInner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference_complement()
            .iter_values()
            .map(render_complement_char(&reference_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(TsArrangementRow::ReferenceComplement, y);
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .reference()
            .iter_values()
            .map(render_source_char(&reference)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(TsArrangementRow::Reference, y);
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query()
            .iter_values()
            .map(render_source_char(&query)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(TsArrangementRow::Query, y);
    y += typewriter::FONT.character_height;

    ts_group = ts_group.add(svg_string(
        ts_arrangement
            .query_complement()
            .iter_values()
            .map(render_complement_char(&query_c)),
        &SvgLocation { x: 0.0, y },
        &typewriter::FONT,
    ));
    rows.insert(TsArrangementRow::QueryComplement, y);
    y += typewriter::FONT.character_height;

    for (identifier, query_inner) in ts_arrangement.query_inners().iter() {
        ts_group = ts_group.add(svg_string(
            query_inner.sequence().iter_values().map(render_inner_char(
                match query_inner.template_switch().primary {
                    TemplateSwitchPrimary::Reference => &reference,
                    TemplateSwitchPrimary::Query => &query,
                },
            )),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        rows.insert(TsArrangementRow::QueryInner { index: identifier }, y);
        y += typewriter::FONT.character_height;
    }

    debug!("Rendering arrows");
    for arrow in &arrows {
        ts_group = ts_group.add(arrow.render(|row| *rows.get(row).unwrap()));
    }

    debug!("Rendering numbers");
    for number in &numbers {
        ts_group = ts_group.add(number.render(|row| *rows.get(row).unwrap()));
    }

    let ts_group_width = ts_arrangement.width() as f32 * typewriter::FONT.character_width;
    let ts_group_height = y;

    svg = svg.add(ts_group);

    let (view_box_width, view_box_height) = if let Some(no_ts_arrangement) = no_ts_arrangement {
        debug!("Rendering no-TS arrangement");

        let vertical_spacer_height = typewriter::FONT.character_height;
        let mut no_ts_group = Group::new().set(
            "transform",
            SvgLocation {
                x: 10.0,
                y: 10.0 + ts_group_height + vertical_spacer_height,
            }
            .as_transform(),
        );
        let mut y = 0.0;

        no_ts_group = no_ts_group.add(svg_string(
            no_ts_arrangement
                .reference()
                .iter_values()
                .map(render_source_char(&reference)),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        y += typewriter::FONT.character_height;

        no_ts_group = no_ts_group.add(svg_string(
            no_ts_arrangement
                .query()
                .iter_values()
                .map(render_source_char(&query)),
            &SvgLocation { x: 0.0, y },
            &typewriter::FONT,
        ));
        y += typewriter::FONT.character_height;

        let no_ts_group_width = no_ts_arrangement.width() as f32 * typewriter::FONT.character_width;
        let no_ts_group_height = y;

        svg = svg.add(no_ts_group);

        (
            ts_group_width.max(no_ts_group_width) + 20.0,
            ts_group_height + vertical_spacer_height + no_ts_group_height + 20.0,
        )
    } else {
        (ts_group_width + 20.0, ts_group_height + 20.0)
    };

    svg = svg.set("viewBox", (0, 0, view_box_width, view_box_height));
    svg::write(output, &svg)?;

    Ok(())
}

fn render_source_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&SourceChar) -> Character<CharacterData> {
    |source_char| match source_char {
        SourceChar::Source {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        SourceChar::Hidden { .. } | SourceChar::Spacer | SourceChar::Blank => {
            Character::new_char_with_default(' ')
        }
    }
}

fn render_complement_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&ComplementChar) -> Character<CharacterData> {
    |complement_char| match complement_char {
        ComplementChar::Complement {
            column,
            lower_case,
            source_hidden,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };

            Character::new_char(
                c,
                if *source_hidden {
                    CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
                } else {
                    Default::default()
                },
            )
        }
        ComplementChar::Gap { source_hidden } => Character::new_char(
            '-',
            if *source_hidden {
                CharacterData::new_colored(COMPLEMENT_SOURCE_HIDDEN_COLOR)
            } else {
                Default::default()
            },
        ),
        ComplementChar::Hidden { .. } | ComplementChar::Blank => {
            Character::new_char(' ', Default::default())
        }
    }
}

fn render_inner_char(
    source_sequence: &IndexedStr,
) -> impl Fn(&InnerChar) -> Character<CharacterData> {
    |inner_char| match inner_char {
        InnerChar::Inner {
            column,
            lower_case,
            copy_depth,
        } => {
            let c = source_sequence.char_at(column.into());
            let c = if *lower_case {
                c.to_ascii_lowercase()
            } else {
                c
            };
            Character::new_char(c, CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Gap { copy_depth } => {
            Character::new_char('-', CharacterData::new_colored(copy_color(copy_depth)))
        }
        InnerChar::Blank => Character::new_char(' ', Default::default()),
    }
}

fn copy_color(copy_depth: &Option<usize>) -> impl ToString {
    if let Some(copy_depth) = copy_depth {
        COPY_COLORS[copy_depth % COPY_COLORS.len()]
    } else {
        "black"
    }
}

pub fn create_error_svg(output: impl Write, error: Error) -> Result<()> {
    let svg = Document::new()
        .set("viewBox", (0, 0, 1920, 1080))
        .add(Circle::new().set("r", 1e5).set("fill", "white"))
        .add(
            Text::new(error.to_string())
                .set("transform", SvgLocation { x: 10.0, y: 10.0 }.as_transform())
                .set("font", "Sans serif"),
        );

    svg::write(output, &svg)?;
    Ok(())
}
