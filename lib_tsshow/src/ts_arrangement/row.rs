use std::fmt::Display;

use super::index_types::TsInnerIdentifier;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TsArrangementRow {
    Reference,
    Query,
    ReferenceComplement,
    QueryComplement,
    Inner { index: TsInnerIdentifier },
}

impl Display for TsArrangementRow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TsArrangementRow::Reference => write!(f, "R"),
            TsArrangementRow::Query => write!(f, "Q"),
            TsArrangementRow::ReferenceComplement => write!(f, "RC"),
            TsArrangementRow::QueryComplement => write!(f, "QC"),
            TsArrangementRow::Inner { index } => write!(f, "I-{index}"),
        }
    }
}
