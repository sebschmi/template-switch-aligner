use std::ops::Range;

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use log::info;

#[derive(Debug, Clone)]
pub struct ChainingSeeds {
    seeds: Vec<ChainingSeed>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ChainingSeed {
    reference_block: Range<usize>,
    query_block: Range<usize>,
}

impl ChainingSeeds {
    /// Compute a set of seeds for the given reference and query sequences.
    ///
    /// The seeds are computed by subdividing the reference sequence into non-overlapping blocks of size `block_size`,
    /// and collecting all their matches in the query sequence.
    /// The last block is merged with the second-to-last block if it is smaller than `block_size`.
    pub fn seed_nonoverlapping<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
        block_size: usize,
    ) -> Self {
        info!("Computing non-overlapping chaining seeds...");
        assert!(
            reference.len() >= block_size,
            "Reference (length: {}) is shorter than the block size {}",
            reference.len(),
            block_size
        );
        assert!(block_size > 0, "Block size must be positive, but is zero");

        // Read into strings so we can use stdlib string matching and don't need to implement our own.
        let reference = reference.as_string();
        let query = query.as_string();

        let reference_block_ranges: Vec<_> =
            nonoverlapping_block_ranges(reference.len(), block_size).collect();
        let reference_blocks: Vec<_> = reference_block_ranges
            .iter()
            .cloned()
            .map(|block_range| &reference[block_range])
            .collect();
        let mut seeds: Vec<_> = find_all_substrings(&query, &reference_blocks)
            .map(
                |SubstringMatch {
                     haystack_offset,
                     needle_index,
                 }| {
                    let reference_block = reference_block_ranges[needle_index].clone();
                    let block_size = reference_block.len();
                    let query_block = haystack_offset..haystack_offset + block_size;

                    ChainingSeed {
                        reference_block,
                        query_block,
                    }
                },
            )
            .collect();
        seeds.sort_unstable();

        ChainingSeeds { seeds }
    }

    pub fn seeds(&self) -> &[ChainingSeed] {
        &self.seeds
    }
}

fn nonoverlapping_block_ranges(
    length: usize,
    block_size: usize,
) -> impl Iterator<Item = Range<usize>> {
    debug_assert!(length >= block_size);
    debug_assert!(block_size > 0);

    (0..length / block_size - 1)
        .map(move |block_index| block_index * block_size..(block_index + 1) * block_size)
        .chain(std::iter::once(
            (length / block_size - 1) * block_size..length,
        ))
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SubstringMatch {
    haystack_offset: usize,
    needle_index: usize,
}

fn find_all_substrings<'haystack, 'needles>(
    haystack: &'haystack str,
    needles: &'needles [&str],
) -> impl use<'haystack, 'needles> + Iterator<Item = SubstringMatch> {
    // Use a slow but simple algorithm.

    struct FindAll<'haystack, 'needle> {
        haystack: &'haystack str,
        haystack_offset: usize,
        needle: &'needle str,
    }

    impl Iterator for FindAll<'_, '_> {
        type Item = usize;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(offset) = self.haystack.find(self.needle) {
                let result = offset + self.haystack_offset;
                self.haystack = &self.haystack[offset + 1..];
                self.haystack_offset += offset + 1;
                Some(result)
            } else {
                let offset = self.haystack.len();
                self.haystack = &self.haystack[offset..];
                self.haystack_offset += offset;
                None
            }
        }
    }

    needles
        .iter()
        .enumerate()
        .flat_map(|(needle_index, needle)| {
            FindAll {
                haystack,
                haystack_offset: 0,
                needle,
            }
            .map(move |haystack_offset| SubstringMatch {
                haystack_offset,
                needle_index,
            })
        })
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

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::{alphabets::dna_alphabet::DnaAlphabet, vec_sequence::VectorGenome},
        interface::sequence::{GenomeSequence, OwnedGenomeSequence},
    };

    use super::{
        find_all_substrings, nonoverlapping_block_ranges, ChainingSeed, ChainingSeeds,
        SubstringMatch,
    };

    #[test]
    fn test_nonoverlapping_block_ranges() {
        assert_eq!(
            nonoverlapping_block_ranges(11, 3).collect::<Vec<_>>(),
            vec![0..3, 3..6, 6..11]
        );
        assert_eq!(
            nonoverlapping_block_ranges(11, 11).collect::<Vec<_>>(),
            vec![0..11]
        );
        assert_eq!(
            nonoverlapping_block_ranges(12, 4).collect::<Vec<_>>(),
            vec![0..4, 4..8, 8..12]
        );
    }

    #[test]
    fn test_find_all_substrings() {
        let haystack = "AAAACATAAA";
        let needles = ["AAA", "ATA", "CAT", "ACA"];
        let mut expected = [(0, 0), (1, 0), (7, 0), (5, 1), (4, 2), (3, 3)].map(
            |(haystack_offset, needle_index)| SubstringMatch {
                haystack_offset,
                needle_index,
            },
        );
        expected.sort();
        let mut actual: Vec<_> = find_all_substrings(haystack, &needles).collect();
        actual.sort();

        assert_eq!(&expected, actual.as_slice());
    }

    #[test]
    fn test_seed_nonoverlapping() {
        let reference = VectorGenome::<DnaAlphabet>::from_slice_u8(b"ACTTGGAAAA").unwrap();
        let query = VectorGenome::<DnaAlphabet>::from_slice_u8(b"TACTGGAAAAACT").unwrap();
        let mut expected = [
            (0..3, 1..4),
            (0..3, 10..13),
            (3..6, 3..6),
            (6..10, 6..10),
            (6..10, 7..11),
        ]
        .map(|(reference_block, query_block)| ChainingSeed {
            reference_block,
            query_block,
        });
        expected.sort();
        let actual = ChainingSeeds::seed_nonoverlapping(
            reference.as_genome_subsequence(),
            query.as_genome_subsequence(),
            3,
        )
        .seeds;

        assert_eq!(&expected, actual.as_slice());
    }
}
