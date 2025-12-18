use std::{cmp::Ordering, iter::Peekable};

use crate::{
    alignment::ts_kind::TsKind,
    anchors::{Anchors, index::AnchorIndex, primary::PrimaryAnchor, secondary::SecondaryAnchor},
};

struct PrimaryAnchorToIndexIter<CoordinateIter: Iterator, AnchorIter: Iterator> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

struct SecondaryAnchorToIndexIter<CoordinateIter: Iterator, AnchorIter: Iterator> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

impl<
    CoordinateIter: Iterator<Item = PrimaryAnchor>,
    AnchorIter: Iterator<Item = (AnchorIndex, PrimaryAnchor)>,
> Iterator for PrimaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Option<AnchorIndex>;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.cmp(anchor) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                    return Some(None);
                }
                Ordering::Equal => {
                    let result = Some(Some(*anchor_index));
                    self.coordinate_iter.next().unwrap();
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        self.coordinate_iter.next().map(|_| None)
    }
}

impl<
    CoordinateIter: Iterator<Item = SecondaryAnchor>,
    AnchorIter: Iterator<Item = (AnchorIndex, SecondaryAnchor)>,
> Iterator for SecondaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Option<AnchorIndex>;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.cmp(anchor) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                    return Some(None);
                }
                Ordering::Equal => {
                    let result = Some(Some(*anchor_index));
                    self.coordinate_iter.next().unwrap();
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        self.coordinate_iter.next().map(|_| None)
    }
}

impl Anchors {
    /// Returns an iterator over the primary anchor indices that correspond to the given primary anchors.
    ///
    /// If a primary anchor does not exist, then the iterator returns `Some(None)`.
    pub fn primary_anchor_to_index_iter(
        &self,
        iter: impl IntoIterator<Item = PrimaryAnchor>,
    ) -> impl Iterator<Item = Option<AnchorIndex>> {
        PrimaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_primaries().peekable(),
        }
    }
    /// Returns an iterator over the secondary anchor indices that correspond to the given secondary anchors.
    ///
    /// If a secondary anchor does not exist, then the iterator returns `Some(None)`.
    pub fn secondary_anchor_to_index_iter(
        &self,
        iter: impl IntoIterator<Item = SecondaryAnchor>,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = Option<AnchorIndex>> {
        SecondaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_secondaries(ts_kind).peekable(),
        }
    }
}
