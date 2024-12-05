use std::ops::Range;

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

#[derive(Debug, Clone)]
pub struct ChainingSeeds {
    seeds: Vec<ChainingSeed>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ChainingSeed {
    reference_block: Range<usize>,
    query_block: Range<usize>,
}

/// Compute a set of seeds for the given reference and query sequences.
///
/// The seeds are computed by subdividing the reference sequence into non-overlapping blocks of size `block_size`,
/// and collecting all their matches in the query sequence.
/// The last block is merged with the second-to-last block if it is smaller than `block_size`.
pub fn seed_nonoverlapping<
    AlphabetType: Alphabet,
    SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType>,
>(
    reference: &SubsequenceType,
    query: &SubsequenceType,
    block_size: usize,
) -> ChainingSeeds {
    assert!(
        reference.len() >= block_size,
        "Reference (length: {}) is shorter than the block size {}",
        reference.len(),
        block_size
    );

    todo!()
}

impl Ord for ChainingSeed {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.reference_block
            .start
            .cmp(&other.reference_block.start)
            .then_with(|| self.reference_block.end.cmp(&other.reference_block.end))
            .then_with(|| self.query_block.start.cmp(&other.query_block.start))
            .then_with(|| self.query_block.end.cmp(&other.query_block.end))
    }
}

impl PartialOrd for ChainingSeed {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
