use super::index_types::SourceColumn;

pub trait Char {
    fn source_column(&self) -> SourceColumn;

    fn is_char(&self) -> bool;

    fn is_gap(&self) -> bool;

    fn is_blank(&self) -> bool;

    fn is_gap_or_blank(&self) -> bool {
        self.is_gap() || self.is_blank()
    }

    fn is_source_char(&self) -> bool;
}
