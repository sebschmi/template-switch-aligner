use std::ops::{Index, IndexMut};

use ndarray::Array2;

use super::BaseAlignmentType;

pub mod iterators;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentMatrixIndex {
    pub(in crate::alignment_matrix) reference_index: usize,
    pub(in crate::alignment_matrix) query_index: usize,
}

impl AlignmentMatrixIndex {
    pub fn new(reference_index: usize, query_index: usize) -> Self {
        Self {
            reference_index,
            query_index,
        }
    }

    pub fn insertion_predecessor(&self) -> Self {
        debug_assert!(self.query_index > 0);

        Self {
            reference_index: self.reference_index,
            query_index: self.query_index - 1,
        }
    }

    pub fn deletion_predecessor(&self) -> Self {
        debug_assert!(self.reference_index > 0);

        Self {
            reference_index: self.reference_index - 1,
            query_index: self.query_index,
        }
    }

    pub fn match_or_substitution_predecessor(&self) -> Self {
        debug_assert!(self.reference_index > 0);
        debug_assert!(self.query_index > 0);

        Self {
            reference_index: self.reference_index - 1,
            query_index: self.query_index - 1,
        }
    }

    pub fn predecessor(&self, alignment_type: BaseAlignmentType) -> Self {
        match alignment_type {
            BaseAlignmentType::None => {
                panic!("Predecessor type 'None' has no predecessor")
            }
            BaseAlignmentType::Insertion => self.insertion_predecessor(),
            BaseAlignmentType::Deletion => self.deletion_predecessor(),
            BaseAlignmentType::Match | BaseAlignmentType::Substitution => {
                self.match_or_substitution_predecessor()
            }
        }
    }
}

impl<T> Index<AlignmentMatrixIndex> for Array2<T> {
    type Output = <Array2<T> as Index<[usize; 2]>>::Output;

    fn index(&self, index: AlignmentMatrixIndex) -> &Self::Output {
        &self[[index.reference_index, index.query_index]]
    }
}

impl<T> IndexMut<AlignmentMatrixIndex> for Array2<T> {
    fn index_mut(&mut self, index: AlignmentMatrixIndex) -> &mut Self::Output {
        &mut self[[index.reference_index, index.query_index]]
    }
}
