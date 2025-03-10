use std::{fmt::Debug, hash::Hash};

use compact_genome::interface::{alphabet::Alphabet, sequence::GenomeSequence};

#[derive(Debug, Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SequencePair {
    reference: String,
    reference_rc: String,
    query: String,
    query_rc: String,
}

impl SequencePair {
    pub fn new<
        AlphabetType: Alphabet,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        reference: &SubsequenceType,
        query: &SubsequenceType,
    ) -> Self {
        Self {
            reference: reference.as_string(),
            reference_rc: reference
                .reverse_complement_iter()
                .map(Into::<char>::into)
                .collect(),
            query: query.as_string(),
            query_rc: query
                .reverse_complement_iter()
                .map(Into::<char>::into)
                .collect(),
        }
    }
}
