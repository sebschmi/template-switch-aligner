use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};
use index::{
    iterators::{
        AlignmentMatrixInnerIterator, AlignmentMatrixQueryIterator,
        AlignmentMatrixReferenceIterator,
    },
    AlignmentMatrixIndex,
};
use ndarray::Array2;

use crate::{alignment_configuration::AlignmentConfiguration, score::Score};

pub mod index;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentMatrix {
    matrix: Array2<AlignmentMatrixEntry>,
    configuration: AlignmentConfiguration,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentMatrixEntry {
    pub score: Score,
    pub alignment_type: BaseAlignmentType,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseAlignmentType {
    /// Marks the matrix origin at [0, 0].
    None,
    Insertion,
    Deletion,
    Match,
    Substitution,
}

impl AlignmentMatrix {
    pub fn new(
        configuration: AlignmentConfiguration,
        reference_length: usize,
        query_length: usize,
    ) -> Self {
        Self {
            matrix: Array2::default((reference_length + 1, query_length + 1)),
            configuration,
        }
    }

    pub fn reference_index_iter(&self, query_index: usize) -> AlignmentMatrixReferenceIterator {
        AlignmentMatrixReferenceIterator::new(query_index, self.matrix.dim().0)
    }

    pub fn query_index_iter(&self, reference_index: usize) -> AlignmentMatrixQueryIterator {
        AlignmentMatrixQueryIterator::new(reference_index, self.matrix.dim().1)
    }

    pub fn inner_index_iter(&self) -> AlignmentMatrixInnerIterator {
        AlignmentMatrixInnerIterator::new(AlignmentMatrixIndex::new(
            self.matrix.dim().0,
            self.matrix.dim().1,
        ))
    }

    pub fn align<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SequenceType,
        query: &SequenceType,
    ) -> Score {
        self.initialise();
        self.align_inner(reference, query);
        self.matrix[[self.matrix.dim().0 - 1, self.matrix.dim().1 - 1]].score
    }

    fn initialise(&mut self) {
        // Initialise matrix origin.
        self.matrix[[0, 0]].score = Score::ZERO;
        self.matrix[[0, 0]].alignment_type = BaseAlignmentType::None;

        // Initialise matrix edges.
        for index in self.reference_index_iter(0).skip(1) {
            self.set_deletion_score(index);
        }
        for index in self.query_index_iter(0).skip(1) {
            self.set_insertion_score(index);
        }
    }

    fn align_inner<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        reference: &SequenceType,
        query: &SequenceType,
    ) {
        for index in self.inner_index_iter() {
            self.set_max_score(index, reference, query);
        }
    }

    fn set_insertion_score(&mut self, index: AlignmentMatrixIndex) {
        self.matrix[index] = self.compute_insertion_entry(index);
    }

    fn set_deletion_score(&mut self, index: AlignmentMatrixIndex) {
        self.matrix[index] = self.compute_deletion_entry(index);
    }

    fn set_max_score<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        index: AlignmentMatrixIndex,
        reference: &SequenceType,
        query: &SequenceType,
    ) {
        // Handle matches and substitutions.
        let mut entry = self.compute_match_or_substitution_entry(index, reference, query);

        // Handle insertions.
        let insertion_entry = self.compute_insertion_entry(index);
        if insertion_entry > entry {
            entry = insertion_entry;
        }

        // Handle deletions.
        let deletion_entry = self.compute_deletion_entry(index);
        if deletion_entry > entry {
            entry = deletion_entry;
        }

        self.matrix[index] = entry;
    }

    fn compute_insertion_entry(&self, index: AlignmentMatrixIndex) -> AlignmentMatrixEntry {
        let alignment_type = BaseAlignmentType::Insertion;
        let predecessor_score = self.matrix[index.predecessor(alignment_type)].score;

        AlignmentMatrixEntry {
            score: predecessor_score + self.configuration.score(alignment_type),
            alignment_type,
        }
    }

    fn compute_deletion_entry(&self, index: AlignmentMatrixIndex) -> AlignmentMatrixEntry {
        let alignment_type = BaseAlignmentType::Deletion;
        let predecessor_score = self.matrix[index.predecessor(alignment_type)].score;

        AlignmentMatrixEntry {
            score: predecessor_score + self.configuration.score(alignment_type),
            alignment_type,
        }
    }

    fn compute_match_or_substitution_entry<
        AlphabetType: Alphabet,
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        index: AlignmentMatrixIndex,
        reference: &SequenceType,
        query: &SequenceType,
    ) -> AlignmentMatrixEntry {
        let alignment_type = if reference[index.reference_index - 1] == query[index.query_index - 1]
        {
            BaseAlignmentType::Match
        } else {
            BaseAlignmentType::Substitution
        };
        let precedessor_score = self.matrix[index.match_or_substitution_predecessor()].score;

        AlignmentMatrixEntry {
            score: precedessor_score + self.configuration.score(alignment_type),
            alignment_type,
        }
    }

    #[cfg(test)]
    fn manual_debug_fill(&mut self, entries: impl IntoIterator<Item = AlignmentMatrixEntry>) {
        let mut entries = entries.into_iter();
        for index in self.inner_index_iter() {
            self.matrix[index] = entries.next().unwrap();
        }
        assert!(entries.next().is_none());
    }
}

impl PartialOrd for AlignmentMatrixEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.score.partial_cmp(&other.score) {
            ord @ (Some(core::cmp::Ordering::Greater | core::cmp::Ordering::Less) | None) => ord,
            ord @ Some(core::cmp::Ordering::Equal) => {
                if self.alignment_type == other.alignment_type {
                    ord
                } else {
                    None
                }
            }
        }
    }
}

impl Default for AlignmentMatrixEntry {
    fn default() -> Self {
        Self {
            score: Score::MIN,
            alignment_type: BaseAlignmentType::None,
        }
    }
}

impl core::fmt::Display for AlignmentMatrix {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut score_column_widths = vec![0; self.matrix.dim().1];
        for reference_index in 0..self.matrix.dim().0 {
            for (query_index, score_column_width) in score_column_widths.iter_mut().enumerate() {
                let score = self.matrix[[reference_index, query_index]].score.as_i64();
                let mut local_score_column_width = 1;

                if score < 0 {
                    local_score_column_width += 1;
                }
                let mut limit = 1;
                for _ in 1..=10 {
                    limit *= 10;
                    if score.abs() >= limit {
                        local_score_column_width += 1;
                    }
                }

                *score_column_width = local_score_column_width.max(*score_column_width);
            }
        }

        for reference_index in 0..self.matrix.dim().0 {
            write!(f, "[ ")?;
            #[allow(clippy::needless_range_loop)]
            for query_index in 0..self.matrix.dim().1 {
                write!(
                    f,
                    "{: >width$}",
                    self.matrix[[reference_index, query_index]].score.as_i64(),
                    width = score_column_widths[query_index],
                )?;
                write!(
                    f,
                    "{} ",
                    match self.matrix[[reference_index, query_index]].alignment_type {
                        BaseAlignmentType::None => "N",
                        BaseAlignmentType::Insertion => "I",
                        BaseAlignmentType::Deletion => "D",
                        BaseAlignmentType::Match => "M",
                        BaseAlignmentType::Substitution => "S",
                    }
                )?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use compact_genome::{
        implementation::bit_vec_sequence_store::BitVectorSequenceStore,
        interface::{alphabet::dna_alphabet::DnaAlphabet, sequence_store::SequenceStore},
    };
    use traitsequence::interface::Sequence;

    use crate::{
        alignment_configuration::AlignmentConfiguration,
        alignment_matrix::{AlignmentMatrixEntry, BaseAlignmentType},
    };

    use super::AlignmentMatrix;

    #[test]
    fn test_simple_alignments() {
        let mut sequence_store = BitVectorSequenceStore::<DnaAlphabet>::new();
        let reference = sequence_store.add_from_slice_u8(b"ACG").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCG").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 1.into());

        let mut manual_matrix = matrix.clone();
        manual_matrix.manual_debug_fill(
            [
                (1, BaseAlignmentType::Match),
                (-1, BaseAlignmentType::Deletion),
                (-3, BaseAlignmentType::Deletion),
                (-1, BaseAlignmentType::Insertion),
                (2, BaseAlignmentType::Match),
                (0, BaseAlignmentType::Deletion),
                (-3, BaseAlignmentType::Insertion),
                (0, BaseAlignmentType::Match),
                (1, BaseAlignmentType::Substitution),
                (-5, BaseAlignmentType::Insertion),
                (-2, BaseAlignmentType::Insertion),
                (1, BaseAlignmentType::Match),
            ]
            .into_iter()
            .map(|(score, alignment_type)| AlignmentMatrixEntry {
                score: score.into(),
                alignment_type,
            }),
        );
        assert_eq!(
            matrix, manual_matrix,
            "matrix:\n{matrix}\nmanual_matrix:\n{manual_matrix}"
        );

        let reference = sequence_store.add_from_slice_u8(b"ACCCGT").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCGT").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 3.into());

        let reference = sequence_store.add_from_slice_u8(b"ACGCCCCCT").unwrap();
        let query = sequence_store.add_from_slice_u8(b"ACCCCCGCT").unwrap();
        let reference = sequence_store.get(&reference);
        let query = sequence_store.get(&query);

        let mut matrix = AlignmentMatrix::new(
            AlignmentConfiguration::default(),
            reference.len(),
            query.len(),
        );
        assert_eq!(matrix.align(reference, query), 5.into());
    }
}
