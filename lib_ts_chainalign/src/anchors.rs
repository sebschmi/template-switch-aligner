use std::fmt::Display;

use lib_tsalign::a_star_aligner::alignment_geometry::AlignmentRange;
use log::trace;

use crate::{
    alignment::{
        coordinates::AlignmentCoordinates,
        sequences::AlignmentSequences,
        ts_kind::{TsAncestor, TsDescendant, TsKind},
    },
    anchors::{
        kmer_matches::find_kmer_matches,
        kmers::{Kmer, KmerStore},
    },
};

pub mod kmer_matches;
pub mod kmers;
#[cfg(test)]
mod tests;

#[derive(Debug, PartialEq, Eq)]
pub struct Anchors {
    pub primary: Vec<PrimaryAnchor>,
    pub secondary_11: Vec<SecondaryAnchor>,
    pub secondary_12: Vec<SecondaryAnchor>,
    pub secondary_21: Vec<SecondaryAnchor>,
    pub secondary_22: Vec<SecondaryAnchor>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct PrimaryAnchor {
    seq1: usize,
    seq2: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub struct SecondaryAnchor {
    /// Ancestor right index in the primary sequence.
    ancestor: usize,
    /// Descendant left index in the primary sequence.
    descendant: usize,
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
        let mut secondary_11: Vec<_> = find_kmer_matches(&s1_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s1.len() - ancestor,
                descendant,
            })
            .collect();
        let mut secondary_12: Vec<_> = find_kmer_matches(&s1_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s1.len() - ancestor,
                descendant,
            })
            .collect();
        let mut secondary_21: Vec<_> = find_kmer_matches(&s2_rc_kmers, &s1_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s2.len() - ancestor,
                descendant,
            })
            .collect();
        let mut secondary_22: Vec<_> = find_kmer_matches(&s2_rc_kmers, &s2_kmers)
            .into_iter()
            .map(|(ancestor, descendant)| SecondaryAnchor {
                ancestor: s2.len() - ancestor,
                descendant,
            })
            .collect();

        // Sort anchors.
        primary.sort_unstable_by_key(|primary_anchor| {
            (
                primary_anchor.seq1.min(primary_anchor.seq2),
                primary_anchor.seq1,
                primary_anchor.seq2,
            )
        });
        secondary_11.sort_unstable_by_key(|secondary_anchor| {
            (
                secondary_anchor.ancestor.min(secondary_anchor.descendant),
                secondary_anchor.ancestor,
                secondary_anchor.descendant,
            )
        });
        secondary_12.sort_unstable_by_key(|secondary_anchor| {
            (
                secondary_anchor.ancestor.min(secondary_anchor.descendant),
                secondary_anchor.ancestor,
                secondary_anchor.descendant,
            )
        });
        secondary_21.sort_unstable_by_key(|secondary_anchor| {
            (
                secondary_anchor.ancestor.min(secondary_anchor.descendant),
                secondary_anchor.ancestor,
                secondary_anchor.descendant,
            )
        });
        secondary_22.sort_unstable_by_key(|secondary_anchor| {
            (
                secondary_anchor.ancestor.min(secondary_anchor.descendant),
                secondary_anchor.ancestor,
                secondary_anchor.descendant,
            )
        });

        Self {
            primary,
            secondary_11,
            secondary_12,
            secondary_21,
            secondary_22,
        }
    }

    pub fn secondary(&self, ts_kind: TsKind) -> &[SecondaryAnchor] {
        match (ts_kind.ancestor, ts_kind.descendant) {
            (TsAncestor::Seq1, TsDescendant::Seq1) => &self.secondary_11,
            (TsAncestor::Seq1, TsDescendant::Seq2) => &self.secondary_12,
            (TsAncestor::Seq2, TsDescendant::Seq1) => &self.secondary_21,
            (TsAncestor::Seq2, TsDescendant::Seq2) => &self.secondary_22,
        }
    }
}

impl PrimaryAnchor {
    pub fn new(seq1: usize, seq2: usize) -> Self {
        Self { seq1, seq2 }
    }

    pub fn start(&self) -> AlignmentCoordinates {
        AlignmentCoordinates::Primary {
            a: self.seq1,
            b: self.seq2,
        }
    }

    pub fn end(&self, k: usize) -> AlignmentCoordinates {
        AlignmentCoordinates::Primary {
            a: self.seq1 + k,
            b: self.seq2 + k,
        }
    }

    pub fn chaining_gaps(&self, second: &Self, k: usize) -> Option<(usize, usize)> {
        let gap_start = self.end(k);
        let gap_end = second.start();
        primary_chaining_gaps(gap_start, gap_end)
    }

    pub fn chaining_gaps_from_start(&self, start: AlignmentCoordinates) -> (usize, usize) {
        let gap_end = self.start();
        primary_chaining_gaps(start, gap_end)
            .unwrap_or_else(|| panic!("self: {self}, start: {start}"))
    }

    pub fn chaining_gaps_to_end(&self, end: AlignmentCoordinates, k: usize) -> (usize, usize) {
        let gap_start = self.end(k);
        primary_chaining_gaps(gap_start, end)
            .unwrap_or_else(|| panic!("self: {self}, end: {end}, k: {k}"))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &SecondaryAnchor,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(k);
        let gap_end = second.start(ts_kind);

        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_start.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => gap_start.primary_ordinate_b().unwrap(),
        };
        let gap_end = gap_end.secondary_ordinate_descendant().unwrap();

        gap_end.checked_sub(gap_start)
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.seq1 + 1 == successor.seq1 && self.seq2 + 1 == successor.seq2
    }
}

fn primary_chaining_gaps(
    gap_start: AlignmentCoordinates,
    gap_end: AlignmentCoordinates,
) -> Option<(usize, usize)> {
    let gap1 = gap_end
        .primary_ordinate_a()
        .unwrap()
        .checked_sub(gap_start.primary_ordinate_a().unwrap())?;
    let gap2 = gap_end
        .primary_ordinate_b()
        .unwrap()
        .checked_sub(gap_start.primary_ordinate_b().unwrap())?;

    Some((gap1, gap2))
}

impl SecondaryAnchor {
    pub fn new(ancestor: usize, descendant: usize) -> Self {
        Self {
            ancestor,
            descendant,
        }
    }

    pub fn start(&self, ts_kind: TsKind) -> AlignmentCoordinates {
        AlignmentCoordinates::Secondary {
            ancestor: self.ancestor,
            descendant: self.descendant,
            ts_kind,
        }
    }

    pub fn end(&self, ts_kind: TsKind, k: usize) -> AlignmentCoordinates {
        AlignmentCoordinates::Secondary {
            ancestor: self.ancestor.checked_sub(k).unwrap(),
            descendant: self.descendant + k,
            ts_kind,
        }
    }

    pub fn chaining_gaps(
        &self,
        second: &Self,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<(usize, usize)> {
        let gap_start = self.end(ts_kind, k);
        let gap_end = second.start(ts_kind);

        let gap1 = gap_start
            .secondary_ordinate_ancestor()
            .unwrap()
            .checked_sub(gap_end.secondary_ordinate_ancestor().unwrap())?;
        let gap2 = gap_end
            .secondary_ordinate_descendant()
            .unwrap()
            .checked_sub(gap_start.secondary_ordinate_descendant().unwrap())?;

        Some((gap1, gap2))
    }

    pub fn chaining_jump_gap(
        &self,
        second: &PrimaryAnchor,
        ts_kind: TsKind,
        k: usize,
    ) -> Option<usize> {
        let gap_start = self.end(ts_kind, k);
        let gap_end = second.start();

        let gap_start = gap_start.secondary_ordinate_descendant().unwrap();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => gap_end.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => gap_end.primary_ordinate_b().unwrap(),
        };

        gap_end.checked_sub(gap_start)
    }

    pub fn chaining_jump_gap_from_start(
        &self,
        start: AlignmentCoordinates,
        ts_kind: TsKind,
    ) -> usize {
        let gap_start = match ts_kind.descendant {
            TsDescendant::Seq1 => start.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => start.primary_ordinate_b().unwrap(),
        };
        let gap_end = self.start(ts_kind).secondary_ordinate_descendant().unwrap();

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn chaining_jump_gap_to_end(
        &self,
        end: AlignmentCoordinates,
        ts_kind: TsKind,
        k: usize,
    ) -> usize {
        let gap_start = self
            .end(ts_kind, k)
            .secondary_ordinate_descendant()
            .unwrap();
        let gap_end = match ts_kind.descendant {
            TsDescendant::Seq1 => end.primary_ordinate_a().unwrap(),
            TsDescendant::Seq2 => end.primary_ordinate_b().unwrap(),
        };

        gap_end.checked_sub(gap_start).unwrap()
    }

    pub fn is_direct_predecessor_of(&self, successor: &Self) -> bool {
        self.ancestor - 1 == successor.ancestor && self.descendant + 1 == successor.descendant
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

        write!(f, "S11: [")?;
        let mut once = true;
        for secondary_anchor in &self.secondary_11 {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{secondary_anchor}")?;
        }
        writeln!(f, "]")?;

        write!(f, "S12: [")?;
        let mut once = true;
        for secondary_anchor in &self.secondary_12 {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{secondary_anchor}")?;
        }
        writeln!(f, "]")?;

        write!(f, "S21: [")?;
        let mut once = true;
        for secondary_anchor in &self.secondary_21 {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{secondary_anchor}")?;
        }
        writeln!(f, "]")?;

        write!(f, "S22: [")?;
        let mut once = true;
        for secondary_anchor in &self.secondary_22 {
            if once {
                once = false;
            } else {
                write!(f, ", ")?;
            }
            write!(f, "{secondary_anchor}")?;
        }
        write!(f, "]")?;

        Ok(())
    }
}

impl Display for PrimaryAnchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.seq1, self.seq2)
    }
}

impl Display for SecondaryAnchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.ancestor, self.descendant)
    }
}

impl From<(usize, usize)> for PrimaryAnchor {
    fn from(value: (usize, usize)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl From<(usize, usize)> for SecondaryAnchor {
    fn from(value: (usize, usize)) -> Self {
        Self::new(value.0, value.1)
    }
}
