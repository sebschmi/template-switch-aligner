use std::{cmp::Ordering, iter::Peekable};

use crate::{
    alignment::ts_kind::TsKind,
    anchors::{Anchors, index::AnchorIndex, primary::PrimaryAnchor, secondary::SecondaryAnchor},
};

struct PrimaryAnchorToIndexIter<
    Anchor,
    CoordinateIter: Iterator<Item = Anchor>,
    AnchorIter: Iterator,
> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

struct SecondaryAnchorToIndexIter<CoordinateIter: Iterator, AnchorIter: Iterator> {
    coordinate_iter: Peekable<CoordinateIter>,
    anchor_iter: Peekable<AnchorIter>,
}

pub trait PartialIntoAnchorIndex {
    type IntoPartSource;
    type IntoTarget;

    fn source_part(&self) -> &Self::IntoPartSource;

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget;
}

impl<
    Anchor: PartialIntoAnchorIndex<IntoPartSource = PrimaryAnchor>,
    CoordinateIter: Iterator<Item = Anchor>,
    AnchorIter: Iterator<Item = (AnchorIndex, PrimaryAnchor)>,
> Iterator for PrimaryAnchorToIndexIter<Anchor, CoordinateIter, AnchorIter>
{
    type Item = Anchor::IntoTarget;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.source_part().cmp(anchor) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                }
                Ordering::Equal => {
                    let result = Some(
                        self.coordinate_iter
                            .next()
                            .unwrap()
                            .partial_into(*anchor_index),
                    );
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        None
    }
}

impl<
    Anchor: PartialIntoAnchorIndex<IntoPartSource = SecondaryAnchor>,
    CoordinateIter: Iterator<Item = Anchor>,
    AnchorIter: Iterator<Item = (AnchorIndex, SecondaryAnchor)>,
> Iterator for SecondaryAnchorToIndexIter<CoordinateIter, AnchorIter>
{
    type Item = Anchor::IntoTarget;

    fn next(&mut self) -> Option<Self::Item> {
        while let (Some(coordinate_anchor), Some((anchor_index, anchor))) =
            (self.coordinate_iter.peek(), self.anchor_iter.peek())
        {
            match coordinate_anchor.source_part().cmp(anchor) {
                Ordering::Less => {
                    self.coordinate_iter.next().unwrap();
                }
                Ordering::Equal => {
                    let result = Some(
                        self.coordinate_iter
                            .next()
                            .unwrap()
                            .partial_into(*anchor_index),
                    );
                    self.anchor_iter.next().unwrap();
                    return result;
                }
                Ordering::Greater => {
                    self.anchor_iter.next().unwrap();
                }
            }
        }

        None
    }
}

impl Anchors {
    /// Returns an iterator over the primary anchor indices that correspond to the given primary anchors.
    ///
    /// If a primary anchor does not exist, then the iterator returns `Some(None)`.
    pub fn primary_anchor_to_index_iter<
        Anchor: PartialIntoAnchorIndex<IntoPartSource = PrimaryAnchor>,
    >(
        &self,
        iter: impl IntoIterator<Item = Anchor>,
    ) -> impl Iterator<Item = Anchor::IntoTarget> {
        PrimaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_primaries().peekable(),
        }
    }

    /// Returns an iterator over the secondary anchor indices that correspond to the given secondary anchors.
    ///
    /// If a secondary anchor does not exist, then the iterator returns `Some(None)`.
    pub fn secondary_anchor_to_index_iter<
        Anchor: PartialIntoAnchorIndex<IntoPartSource = SecondaryAnchor>,
    >(
        &self,
        iter: impl IntoIterator<Item = Anchor>,
        ts_kind: TsKind,
    ) -> impl Iterator<Item = Anchor::IntoTarget> {
        SecondaryAnchorToIndexIter {
            coordinate_iter: iter.into_iter().peekable(),
            anchor_iter: self.enumerate_secondaries(ts_kind).peekable(),
        }
    }
}

impl PartialIntoAnchorIndex for PrimaryAnchor {
    type IntoPartSource = PrimaryAnchor;

    type IntoTarget = AnchorIndex;

    fn source_part(&self) -> &Self::IntoPartSource {
        self
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        target
    }
}

impl PartialIntoAnchorIndex for SecondaryAnchor {
    type IntoPartSource = SecondaryAnchor;

    type IntoTarget = AnchorIndex;

    fn source_part(&self) -> &Self::IntoPartSource {
        self
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        target
    }
}

impl<T> PartialIntoAnchorIndex for (PrimaryAnchor, T) {
    type IntoPartSource = PrimaryAnchor;

    type IntoTarget = (AnchorIndex, T);

    fn source_part(&self) -> &Self::IntoPartSource {
        &self.0
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        (target, self.1)
    }
}

impl<T> PartialIntoAnchorIndex for (SecondaryAnchor, T) {
    type IntoPartSource = SecondaryAnchor;

    type IntoTarget = (AnchorIndex, T);

    fn source_part(&self) -> &Self::IntoPartSource {
        &self.0
    }

    fn partial_into(self, target: AnchorIndex) -> Self::IntoTarget {
        (target, self.1)
    }
}
