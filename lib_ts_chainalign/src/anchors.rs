use std::{fmt::Display, time::Instant};

use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;
use log::{info, trace};

use crate::{
    alignment::{sequences::AlignmentSequences, ts_kind::TsKind},
    anchors::{
        index::AnchorIndex,
        kmer_matches::find_kmer_matches,
        kmers::{Kmer, KmerStore},
        primary::PrimaryAnchor,
        secondary::SecondaryAnchor,
    },
};

pub mod index;
pub mod kmer_matches;
pub mod kmers;
pub mod primary;
pub mod reverse_lookup;
pub mod secondary;
#[cfg(test)]
mod tests;

#[derive(Debug, PartialEq, Eq)]
pub struct Anchors {
    primary: Vec<PrimaryAnchor>,
    secondaries: [Vec<SecondaryAnchor>; 4],
}

impl Anchors {
    pub fn new(
        sequences: &AlignmentSequences,
        range: AlignmentRange,
        k: u32,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self {
        if k <= 8 {
            Self::new_with_kmer_store::<u16>(sequences, range, k, rc_fn)
        } else if k <= 16 {
            Self::new_with_kmer_store::<u32>(sequences, range, k, rc_fn)
        } else if k <= 32 {
            Self::new_with_kmer_store::<u64>(sequences, range, k, rc_fn)
        } else if k <= 64 {
            Self::new_with_kmer_store::<u128>(sequences, range, k, rc_fn)
        } else {
            panic!("Only k-mer sizes up to 64 are supported, but got {k}");
        }
    }

    fn new_with_kmer_store<Store: KmerStore>(
        sequences: &AlignmentSequences,
        range: AlignmentRange,
        k: u32,
        rc_fn: &dyn Fn(u8) -> u8,
    ) -> Self {
        let start_time = Instant::now();

        let k = usize::try_from(k).unwrap();
        let s1 = sequences.seq1();
        let s2 = sequences.seq2();
        let s1_rc: Vec<_> = s1.iter().copied().rev().map(rc_fn).collect();
        let s2_rc: Vec<_> = s2.iter().copied().rev().map(rc_fn).collect();

        let s1_kmer_count =
            (range.reference_limit() - range.reference_offset() + 1).saturating_sub(k);
        let s2_kmer_count = (range.query_limit() - range.query_offset() + 1).saturating_sub(k);

        // Compute k-mers.
        let mut s1_kmers: Vec<_> = (range.reference_offset()
            ..range.reference_offset() + s1_kmer_count)
            .map(|offset| (Kmer::<Store>::from(&s1[offset..offset + k]), offset))
            .collect();
        s1_kmers.sort();
        let s1_kmers = s1_kmers;
        let mut s2_kmers: Vec<_> = (range.query_offset()..range.query_offset() + s2_kmer_count)
            .map(|offset| (Kmer::<Store>::from(&s2[offset..offset + k]), offset))
            .collect();
        s2_kmers.sort();
        let s2_kmers = s2_kmers;
        let mut s1_rc_kmers: Vec<_> = (0..s1_rc.len() + 1 - k)
            .map(|offset| (Kmer::<Store>::from(&s1_rc[offset..offset + k]), offset))
            .collect();
        s1_rc_kmers.sort();
        let s1_rc_kmers = s1_rc_kmers;
        let mut s2_rc_kmers: Vec<_> = (0..s2_rc.len() + 1 - k)
            .map(|offset| (Kmer::<Store>::from(&s2_rc[offset..offset + k]), offset))
            .collect();
        s2_rc_kmers.sort();
        let s2_rc_kmers = s2_rc_kmers;

        trace!("s1_kmers: {s1_kmers:?}");
        trace!("s2_kmers: {s2_kmers:?}");

        // Compute anchors.
        let mut primary: Vec<_> = find_kmer_matches(&s1_kmers, &s2_kmers)
            .into_iter()
            .map(|(seq1, seq2)| PrimaryAnchor { seq1, seq2 })
            .collect();
        let secondary_11: Vec<_> = find_kmer_matches(&s1_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s1.len() - ancestor,
                descendant,
            })
            .collect();
        let secondary_12: Vec<_> = find_kmer_matches(&s1_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s1.len() - ancestor,
                descendant,
            })
            .collect();
        let secondary_21: Vec<_> = find_kmer_matches(&s2_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s2.len() - ancestor,
                descendant,
            })
            .collect();
        let secondary_22: Vec<_> = find_kmer_matches(&s2_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s2.len() - ancestor,
                descendant,
            })
            .collect();
        let mut secondaries = [secondary_11, secondary_12, secondary_21, secondary_22];

        // Sort anchors.
        primary.sort_unstable();
        for secondary in &mut secondaries {
            secondary.sort_unstable();
        }

        let end_time = Instant::now();
        let duration = end_time - start_time;

        info!(
            "Found {} anchors ({} + {} + {} + {} + {}) in {:.0}ms",
            primary.len() + secondaries.iter().map(Vec::len).sum::<usize>(),
            primary.len(),
            secondaries[0].len(),
            secondaries[1].len(),
            secondaries[2].len(),
            secondaries[3].len(),
            duration.as_secs_f64() * 1e3,
        );

        Self {
            primary,
            secondaries,
        }
    }

    pub fn primary(&self, index: AnchorIndex) -> &PrimaryAnchor {
        &self.primary[index.as_usize()]
    }

    pub fn primary_len(&self) -> AnchorIndex {
        self.primary.len().into()
    }

    pub fn enumerate_primaries(&self) -> impl Iterator<Item = (AnchorIndex, PrimaryAnchor)> {
        self.primary
            .iter()
            .copied()
            .enumerate()
            .map(|(index, anchor)| (index.into(), anchor))
    }

    fn secondary_anchor_vec(&self, ts_kind: TsKind) -> &Vec<SecondaryAnchor> {
        &self.secondaries[ts_kind.index()]
    }

    pub fn secondary(&self, index: AnchorIndex, ts_kind: TsKind) -> &SecondaryAnchor {
        &self.secondary_anchor_vec(ts_kind)[index.as_usize()]
    }

    pub fn secondary_len(&self, ts_kind: TsKind) -> AnchorIndex {
        self.secondary_anchor_vec(ts_kind).len().into()
    }

    pub fn enumerate_secondaries(
        &self,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = (AnchorIndex, SecondaryAnchor)> {
        self.secondary_anchor_vec(ts_kind)
            .iter()
            .copied()
            .enumerate()
            .map(|(index, anchor)| (index.into(), anchor))
    }
}

impl Display for Anchors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "P: [")?;
        let mut once = true;
        for primary_anchor in &self.primary {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{primary_anchor}")?;
        }
        writeln!(f, "]")?;

        let mut ts_kind_once = true;
        for ts_kind in TsKind::iter() {
            if ts_kind_once {
                ts_kind_once = false;
            } else {
                writeln!(f)?;
            }

            write!(f, "S{}: [", ts_kind.digits())?;
            let mut once = true;
            for secondary_anchor in self.secondary_anchor_vec(ts_kind) {
                if once {
                    once = false;
                } else {
                    write!(f, ", ")?;
                }
                write!(f, "{secondary_anchor}")?;
            }
            write!(f, "]")?;
        }

        Ok(())
    }
}
