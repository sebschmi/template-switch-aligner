use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("An IO error occurred: {0}.")]
    Io(#[from] std::io::Error),

    #[error("A parsing error occurred on string '{input}': {kind:?}.")]
    Parser {
        input: String,
        kind: nom::error::ErrorKind,
    },

    #[error("Parsing was unsuccessful due to incomplete input: {0:?}.")]
    ParserIncomplete(nom::Needed),

    #[error("The cost table name {0} was encountered twice.")]
    DuplicateCostTableName(String),

    #[error("The template switch cost file contained a wrong set of cost table names. Expected: {expected:?}. Actual: {actual:?}.")]
    WrongCostTableNames {
        actual: Vec<String>,
        expected: Vec<String>,
    },

    #[error("A cost function was attempted to create from a sequence whose index does not strictly increase at {index}.")]
    CostFunctionIndexNotIncreasing { index: usize },
}