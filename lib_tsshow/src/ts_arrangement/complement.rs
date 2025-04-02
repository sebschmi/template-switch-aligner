use super::SourceColumn;

pub struct TsComplementArrangement {
    reference_c: Vec<ComplementChar>,
    query_c: Vec<ComplementChar>,
}

#[derive(Debug, Clone, Copy)]
pub enum ComplementChar {
    Complement {
        column: SourceColumn,
        lower_case: bool,
    },
    Hidden {
        column: SourceColumn,
    },
    Gap,
    Blank,
}

impl TsComplementArrangement {
    pub fn new() -> Self {
        todo!()
    }
}
