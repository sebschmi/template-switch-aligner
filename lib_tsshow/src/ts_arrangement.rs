use character::Char;
use complement::{ComplementChar, TsComplementArrangement};
use index_types::{ArrangementCharColumn, ArrangementColumn, SourceColumn, TsInnerIdentifier};
use inner::{TsInner, TsInnerArrangement};
use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;
use log::debug;
use source::{SourceChar, TsSourceArrangement};
use tagged_vec::TaggedVec;
use template_switch::TemplateSwitch;

use crate::error::Result;

pub mod character;
pub mod complement;
pub mod index_types;
pub mod inner;
pub mod row;
pub mod source;
pub mod template_switch;

pub struct TsArrangement {
    source: TsSourceArrangement,
    complement: TsComplementArrangement,
    inner: TsInnerArrangement,
}

impl TsArrangement {
    pub fn new(
        reference_alignment_offset: usize,
        query_alignment_offset: usize,
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
    ) -> Result<Self> {
        let mut template_switches = Vec::new();
        let mut source = TsSourceArrangement::new(
            reference_alignment_offset,
            query_alignment_offset,
            reference_length,
            query_length,
            alignment,
            &mut template_switches,
        )?;
        let mut complement = TsComplementArrangement::new(&source);
        let inner = TsInnerArrangement::new(&mut source, &mut complement, template_switches);

        Ok(Self {
            source,
            complement,
            inner,
        })
    }

    /// Removes columns that contain only blanks and hidden characters.
    pub fn remove_empty_columns(&mut self) {
        let mut remove_columns = Vec::new();

        'column_iter: for column in self.reference().iter_indices() {
            if !self.reference()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.query()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.reference_complement()[column].is_blank_or_hidden() {
                continue;
            }
            if !self.query_complement()[column].is_blank_or_hidden() {
                continue;
            }
            for inner in self.inner.inners().iter_values() {
                if !inner.sequence()[column].is_blank_or_hidden() {
                    continue 'column_iter;
                }
            }

            remove_columns.push(column);
        }

        let remove_columns = remove_columns;
        debug!("Removing columns: {remove_columns:?}");

        let removed_hidden_chars = self.source.remove_columns(remove_columns.iter().copied());
        self.complement
            .remove_columns(remove_columns.iter().copied());
        self.inner
            .remove_columns(remove_columns.iter().copied(), &removed_hidden_chars);
    }

    /// Unhides all characters of a complement, if any character of the complement is visible.
    pub fn show_complete_complements_if_used(&mut self) {
        let show = |sequence: &mut TaggedVec<ArrangementColumn, ComplementChar>| {
            if sequence.iter_values().any(Char::is_visible_char) {
                sequence.iter_values_mut().for_each(|c| {
                    if c.is_char() {
                        c.show()
                    }
                })
            }
        };

        show(self.complement.reference_complement_mut());
        show(self.complement.query_complement_mut());
    }

    pub fn reference(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        self.source.reference()
    }

    pub fn query(&self) -> &TaggedVec<ArrangementColumn, SourceChar> {
        self.source.query()
    }

    pub fn reference_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        self.complement.reference_complement()
    }

    pub fn query_complement(&self) -> &TaggedVec<ArrangementColumn, ComplementChar> {
        self.complement.query_complement()
    }

    pub fn inners(&self) -> &TaggedVec<TsInnerIdentifier, TsInner> {
        self.inner.inners()
    }

    pub fn reference_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.reference_inners()
    }

    pub fn query_inners(&self) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.query_inners()
    }

    pub fn reference_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.reference_complement_inners()
    }

    pub fn query_complement_inners(
        &self,
    ) -> impl DoubleEndedIterator<Item = (TsInnerIdentifier, &TsInner)> {
        self.inner.query_complement_inners()
    }

    pub fn template_switches(&self) -> impl Iterator<Item = (TsInnerIdentifier, &TemplateSwitch)> {
        self.inners()
            .iter()
            .map(|(identifier, inner)| (identifier, inner.template_switch()))
    }

    pub fn width(&self) -> usize {
        self.source.reference().len()
    }

    pub fn reference_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        self.source
            .reference_arrangement_char_to_arrangement_column(column)
    }

    pub fn query_arrangement_char_to_arrangement_column(
        &self,
        column: ArrangementCharColumn,
    ) -> ArrangementColumn {
        self.source
            .query_arrangement_char_to_arrangement_column(column)
    }

    pub fn reference_source_to_arrangement_column(
        &self,
        column: SourceColumn,
    ) -> ArrangementColumn {
        self.source.reference_source_to_arrangement_column(column)
    }

    pub fn query_source_to_arrangement_column(&self, column: SourceColumn) -> ArrangementColumn {
        self.source.query_source_to_arrangement_column(column)
    }

    pub fn inner_first_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        self.inner.inner_first_non_blank_column(inner_identifier)
    }

    pub fn inner_last_non_blank_column(
        &self,
        inner_identifier: TsInnerIdentifier,
    ) -> ArrangementColumn {
        self.inner.inner_last_non_blank_column(inner_identifier)
    }
}
