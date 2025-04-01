use strong_type::StrongType;

pub struct SourceCoordinates {
    column: SourceColumn,
    row: SourceRow,
}

/// A column index in a source string.
#[derive(StrongType)]
#[strong_type(conversion)]
pub struct SourceColumn(usize);

/// A source row, identifying a source string.
#[derive(StrongType)]
#[strong_type(conversion)]
pub struct SourceRow(usize);

/// A column index in an arrangement.
///
/// This is the index that represents where the character will actually be rendered.
#[derive(StrongType)]
#[strong_type(conversion)]
pub(super) struct ArrangementColumn(usize);

/// A row index in an arrangement.
///
/// This is used to abstract from row keys in the data structure.
#[derive(StrongType)]
#[strong_type(conversion)]
pub(super) struct ArrangementRow(usize);

/// A unique identifier for a gap character.
#[derive(StrongType)]
#[strong_type(conversion)]
pub(super) struct GapIdentifier(usize);

#[derive(Default)]
pub(super) struct GapIdentifierGenerator {
    next: GapIdentifier,
}

impl GapIdentifierGenerator {
    pub fn next(&mut self) -> GapIdentifier {
        let result = self.next;
        self.next = GapIdentifier(self.next.0 + 1);
        result
    }
}

/// The index of an alignment.
#[derive(StrongType)]
#[strong_type(conversion)]
pub(super) struct AlignmentIndex(usize);
