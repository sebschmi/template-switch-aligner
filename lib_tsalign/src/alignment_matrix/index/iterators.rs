use super::AlignmentMatrixIndex;

pub struct AlignmentMatrixReferenceIterator {
    index: AlignmentMatrixIndex,
    limit: usize,
}

pub struct AlignmentMatrixQueryIterator {
    index: AlignmentMatrixIndex,
    limit: usize,
}

/// An iterator over the alignment matrix indices skipping row and column zero.
///
/// The iterator is query-major, i.e. it increments the reference position every iteration, and increments the query position only after reaching the limit of the reference.
pub struct AlignmentMatrixInnerIterator {
    index: AlignmentMatrixIndex,
    limit: AlignmentMatrixIndex,
}

impl AlignmentMatrixReferenceIterator {
    pub(in crate::alignment_matrix) fn new(query_index: usize, limit: usize) -> Self {
        Self {
            index: AlignmentMatrixIndex::new(0, query_index),
            limit,
        }
    }
}

impl AlignmentMatrixQueryIterator {
    pub(in crate::alignment_matrix) fn new(reference_index: usize, limit: usize) -> Self {
        Self {
            index: AlignmentMatrixIndex::new(reference_index, 0),
            limit,
        }
    }
}

impl AlignmentMatrixInnerIterator {
    pub(in crate::alignment_matrix) fn new(limit: AlignmentMatrixIndex) -> Self {
        debug_assert!(limit.reference_index > 0);
        debug_assert!(limit.query_index > 0);

        Self {
            index: if limit.reference_index > 1 && limit.query_index > 1 {
                AlignmentMatrixIndex::new(1, 1)
            } else {
                limit
            },
            limit,
        }
    }
}

impl Iterator for AlignmentMatrixReferenceIterator {
    type Item = AlignmentMatrixIndex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index.reference_index < self.limit {
            let result = self.index;
            self.index.reference_index += 1;
            Some(result)
        } else {
            None
        }
    }
}

impl Iterator for AlignmentMatrixQueryIterator {
    type Item = AlignmentMatrixIndex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index.query_index < self.limit {
            let result = self.index;
            self.index.query_index += 1;
            Some(result)
        } else {
            None
        }
    }
}

impl Iterator for AlignmentMatrixInnerIterator {
    type Item = AlignmentMatrixIndex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index.reference_index < self.limit.reference_index {
            let result = self.index;
            self.index.reference_index += 1;
            Some(result)
        } else if self.index.query_index < self.limit.query_index - 1 {
            self.index.reference_index = 1;
            self.index.query_index += 1;
            let result = self.index;
            self.index.reference_index += 1;
            Some(result)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment_matrix::index::{
        AlignmentMatrixIndex, iterators::AlignmentMatrixInnerIterator,
    };

    #[test]
    fn alignment_matrix_inner_iterator() {
        assert_eq!(
            AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(3, 3)).collect::<Vec<_>>(),
            vec![
                AlignmentMatrixIndex::new(1, 1),
                AlignmentMatrixIndex::new(2, 1),
                AlignmentMatrixIndex::new(1, 2),
                AlignmentMatrixIndex::new(2, 2)
            ],
        );
        assert_eq!(
            AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(1, 3)).collect::<Vec<_>>(),
            vec![],
        );
        assert_eq!(
            AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(3, 1)).collect::<Vec<_>>(),
            vec![],
        );
        assert_eq!(
            AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(1, 1)).collect::<Vec<_>>(),
            vec![],
        );
    }

    #[test]
    #[should_panic]
    fn alignment_matrix_inner_iterator_zero() {
        AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(0, 0));
    }
}
