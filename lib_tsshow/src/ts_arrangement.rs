use character::Char;
use complement::{ComplementChar, TsComplementArrangement};
use index_types::ArrangementColumn;
use inner::{InnerChar, TsInnerArrangement};
use lib_tsalign::a_star_aligner::template_switch_distance::AlignmentType;
use source::{SourceChar, TsSourceArrangement};
use tagged_vec::TaggedVec;
use template_switch::TemplateSwitch;

pub mod character;
pub mod complement;
pub mod index_types;
pub mod inner;
pub mod source;
pub mod template_switch;

pub struct TsArrangement {
    source: TsSourceArrangement,
    complement: TsComplementArrangement,
    inner: TsInnerArrangement,
}

impl TsArrangement {
    pub fn new(
        reference_length: usize,
        query_length: usize,
        alignment: impl IntoIterator<Item = AlignmentType>,
        template_switches_out: &mut Vec<TemplateSwitch>,
    ) -> Self {
        template_switches_out.clear();
        let mut source = TsSourceArrangement::new(
            reference_length,
            query_length,
            alignment,
            template_switches_out,
        );
        let mut complement = TsComplementArrangement::new(&source);
        let inner = TsInnerArrangement::new(&mut source, &mut complement, template_switches_out);

        Self {
            source,
            complement,
            inner,
        }
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
            for inner in self.reference_inners().iter().chain(self.query_inners()) {
                if !inner.1[column].is_blank_or_hidden() {
                    continue 'column_iter;
                }
            }

            remove_columns.push(column);
        }

        self.source.remove_columns(remove_columns.iter().copied());
        self.complement
            .remove_columns(remove_columns.iter().copied());
        self.inner.remove_columns(remove_columns.iter().copied());
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

    pub fn reference_inners(
        &self,
    ) -> &Vec<(TemplateSwitch, TaggedVec<ArrangementColumn, InnerChar>)> {
        self.inner.reference_inners()
    }

    pub fn query_inners(&self) -> &Vec<(TemplateSwitch, TaggedVec<ArrangementColumn, InnerChar>)> {
        self.inner.query_inners()
    }

    pub fn width(&self) -> usize {
        self.source.reference().len()
    }
}
