use super::index_types::SourceColumn;

pub trait Char {
    fn source_column(&self) -> SourceColumn;

    fn is_char(&self) -> bool;

    fn is_gap(&self) -> bool;

    fn is_spacer(&self) -> bool;

    fn is_blank(&self) -> bool;

    fn is_gap_or_blank(&self) -> bool {
        self.is_gap() || self.is_blank()
    }

    fn is_source_char(&self) -> bool;

    /// True if is a character (i.e. no gap) and is visible.
    fn is_visible_char(&self) -> bool {
        self.is_char() && !self.is_hidden()
    }

    fn is_hidden(&self) -> bool;

    fn is_blank_or_hidden(&self) -> bool {
        self.is_blank() || self.is_hidden()
    }

    fn is_blank_or_hidden_or_spacer(&self) -> bool {
        self.is_blank() || self.is_hidden() || self.is_spacer()
    }
}
