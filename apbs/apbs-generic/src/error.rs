// APBS error types
// All C return codes become Result<T, ApbsError>

use thiserror::Error;

#[derive(Error, Debug, Clone)]
pub enum ApbsError {
    #[error("I/O error: {0}")]
    Io(String),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Unsupported format: {0}")]
    UnsupportedFormat(String),

    #[error("Grid error: {0}")]
    Grid(String),

    #[error("Solver error: {0}")]
    Solver(String),

    #[error("Molecule error: {0}")]
    Molecule(String),

    #[error("Memory error: {0}")]
    Memory(String),

    #[error("Assertion failed: {0}")]
    Assertion(String),

    #[error("Input file error: {0}")]
    InputFile(String),
}

pub type ApbsResult<T> = Result<T, ApbsError>;
