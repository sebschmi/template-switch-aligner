use std::{
    io::{Read, Write},
    mem,
    ops::{Index, IndexMut},
    slice,
};

use bitvec::{bitvec, order::LocalBits, vec::BitVec};

pub struct CostArray2D<Cost> {
    len: [usize; 2],
    data: Vec<Cost>,
    is_exact: BitVec<usize, LocalBits>,
}

impl<Cost> CostArray2D<Cost> {
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

    #[expect(dead_code)]
    pub fn write(&self, mut write: impl Write) -> std::io::Result<()>
    where
        Cost: Copy,
    {
        write.write_all(&self.len[0].to_ne_bytes())?;
        write.write_all(&self.len[1].to_ne_bytes())?;

        let cost_size = mem::size_of::<Cost>();
        let data: &[u8] = unsafe {
            slice::from_raw_parts(self.data.as_ptr() as *const u8, self.data.len() * cost_size)
        };
        write.write_all(data)
    }

    #[expect(dead_code)]
    pub fn read(mut read: impl Read) -> std::io::Result<Self>
    where
        Cost: Copy,
    {
        let mut buffer = [0; mem::size_of::<usize>()];
        read.read_exact(&mut buffer)?;
        let len1 = usize::from_ne_bytes(buffer);
        read.read_exact(&mut buffer)?;
        let len2 = usize::from_ne_bytes(buffer);
        let len = [len1, len2];

        let cost_size = mem::size_of::<Cost>();
        let data_len_bytes = cost_size * len1 * len2;

        let mut data = Vec::<Cost>::with_capacity(data_len_bytes);
        let mut data_bytes = unsafe {
            Vec::from_raw_parts(data.as_mut_ptr() as *mut u8, 0, data.capacity() * cost_size)
        };
        read.by_ref()
            .take(data_len_bytes.try_into().unwrap())
            .read_to_end(&mut data_bytes)?;
        unsafe {
            data.set_len(len1 * len2);
        };
        data_bytes.leak();

        Ok(Self {
            len,
            data,
            is_exact: bitvec![usize, LocalBits; 0; len[0] * len[1]],
        })
    }
}

impl<Cost> Index<[usize; 2]> for CostArray2D<Cost> {
    type Output = Cost;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

impl<Cost> IndexMut<[usize; 2]> for CostArray2D<Cost> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.data[coordinates_to_index(index[0], index[1], self.len)]
    }
}

#[inline]
fn coordinates_to_index(c1: usize, c2: usize, len: [usize; 2]) -> usize {
    c1 * len[1] + c2
}
