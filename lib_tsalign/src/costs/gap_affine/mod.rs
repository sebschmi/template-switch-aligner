use std::marker::PhantomData;

use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};

use crate::costs::cost::Cost;

pub mod io;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct GapAffineAlignmentCostTable<AlphabetType> {
    name: String,
    substitution_cost_table: Vec<Cost>,
    gap_open_cost_vector: Vec<Cost>,
    gap_extend_cost_vector: Vec<Cost>,
    phantom_data: PhantomData<AlphabetType>,
}

impl<AlphabetType: Alphabet> GapAffineAlignmentCostTable<AlphabetType> {
    pub fn new(
        name: impl Into<String>,
        substitution_cost_table: impl Into<Vec<Cost>>,
        gap_open_cost_vector: impl Into<Vec<Cost>>,
        gap_extend_cost_vector: impl Into<Vec<Cost>>,
    ) -> Self {
        Self {
            name: name.into(),
            substitution_cost_table: substitution_cost_table.into(),
            gap_open_cost_vector: gap_open_cost_vector.into(),
            gap_extend_cost_vector: gap_extend_cost_vector.into(),
            phantom_data: Default::default(),
        }
    }

    pub fn new_zero() -> Self {
        Self {
            name: "new_zero".to_string(),
            substitution_cost_table: vec![Cost::ZERO; AlphabetType::SIZE * AlphabetType::SIZE],
            gap_open_cost_vector: vec![Cost::ZERO; AlphabetType::SIZE],
            gap_extend_cost_vector: vec![Cost::ZERO; AlphabetType::SIZE],
            phantom_data: Default::default(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn match_or_substitution_cost(
        &self,
        c1: impl Into<AlphabetType::CharacterType>,
        c2: impl Into<AlphabetType::CharacterType>,
    ) -> Cost {
        let c1 = c1.into().index();
        let c2 = c2.into().index();

        self.substitution_cost_table[c1 * AlphabetType::SIZE + c2]
    }

    pub fn gap_open_cost(&self, c: impl Into<AlphabetType::CharacterType>) -> Cost {
        self.gap_open_cost_vector[c.into().index()]
    }

    pub fn gap_extend_cost(&self, c: impl Into<AlphabetType::CharacterType>) -> Cost {
        self.gap_extend_cost_vector[c.into().index()]
    }

    pub fn gap_costs(&self, c: impl Into<AlphabetType::CharacterType>, is_first: bool) -> Cost {
        if is_first {
            self.gap_open_cost(c)
        } else {
            self.gap_extend_cost(c)
        }
    }
}
