use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("An IO error occurred: {0}")]
    Io(#[from] std::io::Error),

    #[error("A parsing error occurred on string '{input}': {kind:?}")]
    Parser {
        input: String,
        kind: nom::error::ErrorKind,
    },

    #[error("Parsing was unsuccessful due to incomplete input: {0:?}")]
    ParserIncomplete(nom::Needed),
}
