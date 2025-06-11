use std::{iter, marker::PhantomData};

use compact_genome::interface::alphabet::{Alphabet, AlphabetCharacter};
use generic_a_star::cost::AStarCost;

pub mod io;

#[derive(Debug, Eq, PartialEq)]
pub struct GapAffineAlignmentCostTable<AlphabetType, Cost> {
    name: String,
    substitution_cost_table: Vec<Cost>,
    gap_open_cost_vector: Vec<Cost>,
    gap_extend_cost_vector: Vec<Cost>,
    phantom_data: PhantomData<AlphabetType>,

    // Cache some values that are potentially used often.
    min_substitution_cost: Cost,
    min_gap_open_cost: Cost,
}

impl<AlphabetType: Alphabet, Cost: AStarCost> GapAffineAlignmentCostTable<AlphabetType, Cost> {
    pub fn new(
        name: impl Into<String>,
        substitution_cost_table: impl Into<Vec<Cost>>,
        gap_open_cost_vector: impl Into<Vec<Cost>>,
        gap_extend_cost_vector: impl Into<Vec<Cost>>,
    ) -> Self {
        let mut result = Self {
            name: name.into(),
            substitution_cost_table: substitution_cost_table.into(),
            gap_open_cost_vector: gap_open_cost_vector.into(),
            gap_extend_cost_vector: gap_extend_cost_vector.into(),
            phantom_data: Default::default(),

            min_substitution_cost: Cost::max_value(),
            min_gap_open_cost: Cost::max_value(),
        };

        result.min_substitution_cost = AlphabetType::iter()
            .flat_map(|c1| {
                iter::repeat(c1)
                    .zip(AlphabetType::iter())
                    .filter_map(|(c1, c2)| {
                        if c1 != c2 {
                            Some(result.substitution_cost(c1.clone(), c2))
                        } else {
                            None
                        }
                    })
            })
            .min()
            .unwrap();
        result.min_gap_open_cost = result.gap_open_cost_vector.iter().min().copied().unwrap();

        result
    }

    pub fn new_base_agnostic(
        name: impl Into<String>,
        match_cost: Cost,
        substitution_cost: Cost,
        gap_open_cost: Cost,
        gap_extend_cost: Cost,
    ) -> Self {
        Self::new(
            name.into(),
            (0..AlphabetType::SIZE)
                .flat_map(|a| {
                    (0..AlphabetType::SIZE).map(move |b| {
                        if a == b {
                            match_cost
                        } else {
                            substitution_cost
                        }
                    })
                })
                .collect::<Vec<_>>(),
            iter::repeat_n(gap_open_cost, AlphabetType::SIZE.into()).collect::<Vec<_>>(),
            iter::repeat_n(gap_extend_cost, AlphabetType::SIZE.into()).collect::<Vec<_>>(),
        )
    }

    pub fn new_zero() -> Self {
        let alphabet_size: usize = AlphabetType::SIZE.into();

        Self::new(
            "new_zero".to_string(),
            vec![Cost::zero(); alphabet_size * alphabet_size],
            vec![Cost::zero(); alphabet_size],
            vec![Cost::zero(); alphabet_size],
        )
    }

    /// Creates a new table with all costs set to `Cost::max_value()`.
    pub fn new_max() -> Self {
        let alphabet_size: usize = AlphabetType::SIZE.into();

        Self::new(
            "new_max".to_string(),
            vec![Cost::max_value(); alphabet_size * alphabet_size],
            vec![Cost::max_value(); alphabet_size],
            vec![Cost::max_value(); alphabet_size],
        )
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
        let c1: usize = c1.into().index().into();
        let c2: usize = c2.into().index().into();

        self.substitution_cost_table[c1 * usize::from(AlphabetType::SIZE) + c2]
    }

    pub fn min_match_cost(&self) -> Cost {
        AlphabetType::iter()
            .map(|character| self.match_cost(character.clone(), character))
            .min()
            .unwrap()
    }

    pub fn min_substitution_cost(&self) -> Cost {
        self.min_substitution_cost
    }

    pub fn gap_open_cost(&self, c: impl Into<AlphabetType::CharacterType>) -> Cost {
        self.gap_open_cost_vector[usize::from(c.into().index())]
    }

    pub fn gap_extend_cost(&self, c: impl Into<AlphabetType::CharacterType>) -> Cost {
        self.gap_extend_cost_vector[usize::from(c.into().index())]
    }

    pub fn gap_costs(&self, c: impl Into<AlphabetType::CharacterType>, is_first: bool) -> Cost {
        if is_first {
            self.gap_open_cost(c)
        } else {
            self.gap_extend_cost(c)
        }
    }

    pub fn min_gap_open_cost(&self) -> Cost {
        self.min_gap_open_cost
    }

    pub fn max_gap_open_cost(&self) -> Cost {
        self.gap_open_cost_vector.iter().max().copied().unwrap()
    }

    pub fn min_gap_extend_cost(&self) -> Cost {
        self.gap_extend_cost_vector.iter().min().copied().unwrap()
    }

    /// The minimum cost caused by a non-matching alignment.
    ///
    /// This assumes that the alignment is not preceded by a gap,
    /// i.e. a gap in the alignment itself would be priced with gap-open costs.
    pub fn min_non_match_cost(&self) -> Cost {
        self.min_gap_open_cost().min(self.min_substitution_cost())
    }

    /// Fill all costs with their minimum over all characters.
    ///
    /// Gap open costs and gap extend costs are set to the minimum value over all characters.
    /// Match costs are set to the minimum value over all matches and substitution costs are set to the minimum value over all substitutions.
    pub fn into_lower_bound(self) -> Self {
        let min_match_cost = self.min_match_cost();
        let min_substitution_cost = self.min_substitution_cost();
        let substitution_cost_table: Vec<_> = AlphabetType::iter()
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

        Self::new(
            self.name,
            substitution_cost_table,
            vec_into_min(self.gap_open_cost_vector),
            vec_into_min(self.gap_extend_cost_vector),
        )
    }

    /// Fill all costs with their minimum over all characters.
    ///
    /// Gap open costs and gap extend costs are set to the minimum value over all characters.
    /// Match and substitution costs are set to the minimum value over all matches and substitutions.
    pub fn into_match_agnostic_lower_bound(self) -> Self {
        Self::new(
            self.name,
            vec_into_min(self.substitution_cost_table),
            vec_into_min(self.gap_open_cost_vector),
            vec_into_min(self.gap_extend_cost_vector),
        )
    }
}

fn vec_into_min<ValueType: Clone + Ord>(mut vec: Vec<ValueType>) -> Vec<ValueType> {
    let min = vec.iter().min().unwrap().clone();
    vec.iter_mut().for_each(|value| *value = min.clone());
    vec
}

impl<AlphabetType, Cost: Clone> Clone for GapAffineAlignmentCostTable<AlphabetType, Cost> {
    fn clone(&self) -> Self {
        Self {
            name: self.name.clone(),
            substitution_cost_table: self.substitution_cost_table.clone(),
            gap_open_cost_vector: self.gap_open_cost_vector.clone(),
            gap_extend_cost_vector: self.gap_extend_cost_vector.clone(),
            phantom_data: self.phantom_data,

            min_substitution_cost: self.min_substitution_cost.clone(),
            min_gap_open_cost: self.min_gap_open_cost.clone(),
        }
    }
}
