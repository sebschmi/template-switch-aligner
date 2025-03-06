use crate::alignment_matrix::BaseAlignmentType;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentConfiguration<Cost> {
    pub match_cost: Cost,
    pub substitution_cost: Cost,
    pub insertion_cost: Cost,
    pub deletion_cost: Cost,
}

impl<Cost: Copy> AlignmentConfiguration<Cost> {
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

impl<Cost: From<u8>> Default for AlignmentConfiguration<Cost> {
    fn default() -> Self {
        Self {
            match_cost: 0.into(),
            substitution_cost: 2.into(),
            insertion_cost: 3.into(),
            deletion_cost: 3.into(),
        }
    }
}
