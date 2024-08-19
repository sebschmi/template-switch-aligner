use crate::{alignment_matrix::BaseAlignmentType, score::Score};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentConfiguration {
    match_score: Score,
    substitution_score: Score,
    insertion_score: Score,
    deletion_score: Score,
}

impl AlignmentConfiguration {
    pub fn score(&self, alignment_type: BaseAlignmentType) -> Score {
        match alignment_type {
            BaseAlignmentType::None => {
                panic!("Alignment type 'None' has no score")
            }
            BaseAlignmentType::Insertion => self.insertion_score,
            BaseAlignmentType::Deletion => self.deletion_score,
            BaseAlignmentType::Match => self.match_score,
            BaseAlignmentType::Substitution => self.substitution_score,
        }
    }
}

impl Default for AlignmentConfiguration {
    fn default() -> Self {
        Self {
            match_score: 1.into(),
            substitution_score: (-1).into(),
            insertion_score: (-2).into(),
            deletion_score: (-2).into(),
        }
    }
}
