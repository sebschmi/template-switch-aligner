use crate::{alignment_matrix::BaseAlignmentType, cost::Cost};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentConfiguration {
    pub match_cost: Cost,
    pub substitution_cost: Cost,
    pub insertion_cost: Cost,
    pub deletion_cost: Cost,
}

impl AlignmentConfiguration {
    pub fn cost(&self, alignment_type: BaseAlignmentType) -> Cost {
        match alignment_type {
            BaseAlignmentType::None => {
                panic!("Alignment type 'None' has no cost")
            }
            BaseAlignmentType::Insertion => self.insertion_cost,
            BaseAlignmentType::Deletion => self.deletion_cost,
            BaseAlignmentType::Match => self.match_cost,
            BaseAlignmentType::Substitution => self.substitution_cost,
        }
    }
}

impl Default for AlignmentConfiguration {
    fn default() -> Self {
        Self {
            match_cost: 0.into(),
            substitution_cost: 2.into(),
            insertion_cost: 3.into(),
            deletion_cost: 3.into(),
        }
    }
}
