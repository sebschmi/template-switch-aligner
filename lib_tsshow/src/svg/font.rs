use std::{borrow::Borrow, collections::HashMap, fmt::Debug};

use svg::node::element::{Group, Path};

use crate::plain_text::mutlipair_alignment_renderer::{Character, NoCharacterData};

use super::SvgLocation;

pub mod sans_serif_mono;
pub mod typewriter;

#[derive(Debug)]
pub struct Font {
    characters: HashMap<char, String>,
    pub character_width: f32,
    pub character_height: f32,
    stroke_width: f32,
    scale: f32,
}

#[derive(Debug, Clone)]
pub struct CharacterData {
    pub color: String,
}

impl Font {
    pub fn new(
        characters: HashMap<char, String>,
        character_width: f32,
        character_height: f32,
        stroke_width: f32,
        scale: f32,
    ) -> Self {
        Self {
            characters,
            character_width: character_width * scale,
            character_height: character_height * scale,
            stroke_width,
            scale,
        }
    }
}

pub fn svg_character<Data: Debug>(
    character: &Character<Data>,
    location: &SvgLocation,
    font: &Font,
) -> Path
where
    CharacterData: for<'a> From<&'a Data>,
{
    let character_path = font
        .characters
        .get(&character.as_char())
        .unwrap_or_else(|| panic!("Unsupported character: {character:?}"));
    let data: CharacterData = character.data().into();
    Path::new()
        .set("fill", data.color.clone())
        .set("stroke", data.color.clone())
        .set("stroke-width", font.stroke_width)
        .set(
            "transform",
            format!("{} scale({})", location.as_transform(), font.scale),
        )
        .set("d", character_path.as_str())
}

pub fn svg_string<
    'character,
    Data: 'character + Debug,
    Item: 'character + Borrow<Character<Data>>,
>(
    string: impl IntoIterator<Item = Item>,
    location: &SvgLocation,
    font: &Font,
) -> Group
where
    CharacterData: for<'a> From<&'a Data>,
{
    let mut group = Group::new();
    group = group.set("transform", location.as_transform());

    for (index, c) in string.into_iter().enumerate() {
        let c = c.borrow();
        group = group.add(svg_character(
            c,
            &SvgLocation {
                x: index as f32 * font.character_width,
                y: 0.0,
            },
            font,
        ));
    }
    group
}

impl CharacterData {
    pub fn new_colored(color: impl ToString) -> Self {
        #[allow(clippy::needless_update)]
        Self {
            color: color.to_string(),
            ..Default::default()
        }
    }
}

impl Default for CharacterData {
    fn default() -> Self {
        Self {
            color: "black".to_string(),
        }
    }
}

impl<'a> From<&'a NoCharacterData> for CharacterData {
    fn from(_: &'a NoCharacterData) -> Self {
        Self::default()
    }
}

impl<'a> From<&'a CharacterData> for CharacterData {
    fn from(value: &'a CharacterData) -> Self {
        value.clone()
    }
}
