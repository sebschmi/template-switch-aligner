use std::ops::{Index, IndexMut};

use bitvec::{bitvec, order::LocalBits, vec::BitVec};

pub struct ChainingCostArray<Cost> {
    len: [usize; 2],
    data: Vec<Cost>,
    is_exact: BitVec<usize, LocalBits>,
}

impl<Cost> ChainingCostArray<Cost> {
    pub fn new_from_cost(len: [usize; 2], cost: Cost) -> Self
    where
        Cost: Clone,
    {
        Self {
            len,
            data: vec![cost; len[0] * len[1]],
            is_exact: bitvec![usize, LocalBits; 0; len[0] * len[1]],
        }
    }

    pub fn is_exact(&self, c1: usize, c2: usize) -> bool {
        self.is_exact[coordinates_to_index(c1, c2, self.len)]
    }

    pub fn set_exact(&mut self, c1: usize, c2: usize) {
        debug_assert!(!self.is_exact(c1, c2));
        self.is_exact
            .set(coordinates_to_index(c1, c2, self.len), true);
    }

    pub fn dim(&self) -> (usize, usize) {
        (self.len[0], self.len[1])
    }
}

impl<Cost> Index<[usize; 2]> for ChainingCostArray<Cost> {
    type Output = Cost;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

impl<Cost> IndexMut<[usize; 2]> for ChainingCostArray<Cost> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

#[inline]
fn coordinates_to_index(c1: usize, c2: usize, len: [usize; 2]) -> usize {
    c1 * len[1] + c2
}
