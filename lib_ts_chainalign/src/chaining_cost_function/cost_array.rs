use std::ops::{Index, IndexMut};

use bitvec::{bitvec, order::LocalBits, vec::BitVec};

pub struct ChainingCostArray<Cost> {
    len: [usize; 2],
    data: Vec<Cost>,
    /// For each ordinate 0, store the indices of the corresponding ordinates 1, ordered by their cost.
    cost_order_permutation: Vec<Vec<u32>>,
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
            cost_order_permutation: vec![Vec::new(); len[0]],
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

    /// Iterate over all ordinates `c2` that belong to this `c1` in order of their cost.
    pub fn iter_in_cost_order(&mut self, c1: usize) -> impl Iterator<Item = (usize, Cost)>
    where
        Cost: Copy + Ord,
    {
        if self.cost_order_permutation[c1].is_empty() {
            self.cost_order_permutation[c1].extend((0..).take(self.len[1]));
            self.cost_order_permutation[c1].sort_unstable_by_key(|c2| {
                self.data[coordinates_to_index(c1, usize::try_from(*c2).unwrap(), self.len)]
            });
        }

        let data = &self.data;
        let len = self.len;
        self.cost_order_permutation[c1]
            .iter()
            .copied()
            .filter_map(move |c2| {
                let c2 = usize::try_from(c2).unwrap();
                (c1 != c2).then(|| (c2, data[coordinates_to_index(c1, c2, len)]))
            })
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
        self.cost_order_permutation[index[0]].clear();
        &mut self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

#[inline]
fn coordinates_to_index(c1: usize, c2: usize, len: [usize; 2]) -> usize {
    c1 * len[1] + c2
}
