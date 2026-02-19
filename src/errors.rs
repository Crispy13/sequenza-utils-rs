use thiserror::Error;

pub type Result<T> = std::result::Result<T, AppError>;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("missing value for argument: {flag}")]
    MissingValue { flag: String },
    #[error("missing required argument: {field}")]
    MissingRequired { field: String },
    #[error("invalid value for {flag}={value}: {reason}")]
    InvalidValue {
        flag: String,
        value: String,
        reason: String,
    },
    #[error("unsupported argument: {arg}")]
    UnsupportedArgument { arg: String },
    #[error("required command not found in PATH: {command}")]
    CommandNotFound { command: String },
    #[error("command failed: {command} (exit: {code:?}) stderr: {stderr}")]
    CommandFailed {
        command: String,
        code: Option<i32>,
        stderr: String,
    },
    #[error("parse error: {message}")]
    ParseError { message: String },
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}
