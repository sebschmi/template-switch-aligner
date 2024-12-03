use std::marker::PhantomData;

use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};

use crate::costs::cost::Cost;

pub mod io;

#[derive(Debug, Eq, PartialEq)]
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

    /// Creates a new table with all costs set to `Cost::MAX`.
    pub fn new_max() -> Self {
        Self {
            name: "new_max".to_string(),
            substitution_cost_table: vec![Cost::MAX; AlphabetType::SIZE * AlphabetType::SIZE],
            gap_open_cost_vector: vec![Cost::MAX; AlphabetType::SIZE],
            gap_extend_cost_vector: vec![Cost::MAX; AlphabetType::SIZE],
            phantom_data: Default::default(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn match_cost(
        &self,
        c1: impl Into<AlphabetType::CharacterType>,
        c2: impl Into<AlphabetType::CharacterType>,
    ) -> Cost {
        let c1 = c1.into();
        let c2 = c2.into();
        debug_assert_eq!(c1.index(), c2.index());

        self.match_or_substitution_cost(c1, c2)
    }

    pub fn substitution_cost(
        &self,
        c1: impl Into<AlphabetType::CharacterType>,
        c2: impl Into<AlphabetType::CharacterType>,
    ) -> Cost {
        let c1 = c1.into();
        let c2 = c2.into();
        debug_assert_ne!(c1.index(), c2.index());

        self.match_or_substitution_cost(c1, c2)
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

    pub fn min_match_cost(&self) -> Cost {
        AlphabetType::iter()
            .map(|character| self.match_cost(character.clone(), character))
            .min()
            .unwrap()
    }

    pub fn min_substitution_cost(&self) -> Cost {
        AlphabetType::iter()
            .flat_map(|c1| {
                AlphabetType::iter().filter_map(move |c2| {
                    if c1 != c2 {
                        Some(self.substitution_cost(c1.clone(), c2))
                    } else {
                        None
                    }
                })
            })
            .min()
            .unwrap()
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

    pub fn min_gap_open_cost(&self) -> Cost {
        self.gap_open_cost_vector.iter().min().copied().unwrap()
    }

    pub fn min_gap_extend_cost(&self) -> Cost {
        self.gap_extend_cost_vector.iter().min().copied().unwrap()
    }

    /// Fill all costs with their minimum over all characters.
    ///
    /// Gap open costs and gap extend costs are set to the minimum value over all characters.
    /// Match costs are set to the minimum value over all matches and substitution costs are set to the minimum value over all substitutions.
    pub fn into_lower_bound(self) -> Self {
        let min_match_cost = self.min_match_cost();
        let min_substitution_cost = self.min_substitution_cost();
        let substitution_cost_table = AlphabetType::iter()
            .flat_map(|c1| {
                AlphabetType::iter().map(move |c2| {
                    if c1 == c2 {
                        min_match_cost
                    } else {
                        min_substitution_cost
                    }
                })
            })
            .collect();

        Self {
            name: self.name,
            substitution_cost_table,
            gap_open_cost_vector: vec_into_min(self.gap_open_cost_vector),
            gap_extend_cost_vector: vec_into_min(self.gap_extend_cost_vector),
            phantom_data: self.phantom_data,
        }
    }

    /// Fill all costs with their minimum over all characters.
    ///
    /// Gap open costs and gap extend costs are set to the minimum value over all characters.
    /// Match and substitution costs are set to the minimum value over all matches and substitutions.
    pub fn into_match_agnostic_lower_bound(self) -> Self {
        Self {
            name: self.name,
            substitution_cost_table: vec_into_min(self.substitution_cost_table),
            gap_open_cost_vector: vec_into_min(self.gap_open_cost_vector),
            gap_extend_cost_vector: vec_into_min(self.gap_extend_cost_vector),
            phantom_data: self.phantom_data,
        }
    }
}

fn vec_into_min<ValueType: Clone + Ord>(mut vec: Vec<ValueType>) -> Vec<ValueType> {
    let min = vec.iter().min().unwrap().clone();
    vec.iter_mut().for_each(|value| *value = min.clone());
    vec
}

impl<AlphabetType> Clone for GapAffineAlignmentCostTable<AlphabetType> {
    fn clone(&self) -> Self {
        Self {
            name: self.name.clone(),
            substitution_cost_table: self.substitution_cost_table.clone(),
            gap_open_cost_vector: self.gap_open_cost_vector.clone(),
            gap_extend_cost_vector: self.gap_extend_cost_vector.clone(),
            phantom_data: self.phantom_data,
        }
    }
}
