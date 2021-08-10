//! Define error type

/// Error enum
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Failled to open fastq file
    #[error("Can't open fastq file {path}")]
    FastqOpenError { path: String },

    /// Fastq parsing error
    #[error("Error durring fastq record parsing")]
    FastqParsingError,

    /// Failled to open bam file
    #[error("Can't open bam file {path}")]
    BamOpenError { path: String },

    /// Bam parsing error
    #[error("Error durring bam record parsing")]
    BamParsingError,
}
