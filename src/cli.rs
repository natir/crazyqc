//! All stuff relate to command line

/// Struct use to parse command line argument
#[derive(clap::Clap, Debug)]
#[clap(
    name = "crazyqc",
    version = "0.1 ",
    author = "Pierre Marijon <pierre.marijon-ext@aphp.fr>",
    about = "A crazily fast tools to compute small set of sequence quality control"
)]
pub struct Command {
    /// Paths of the fastq files
    #[clap(short = 'q', long = "fastq", about = "Fastq input")]
    pub fastq: Vec<String>,

    /// Paths of the bam files
    #[clap(short = 'b', long = "bam", about = "Bam input, optional")]
    pub bam: Option<Vec<String>>,

    /// Paths where result will be written
    #[clap(
        short = 'o',
        long = "output",
        about = "Path where result will be write, default: stdout"
    )]
    pub output: Option<String>,

    /// Control read buffer size
    #[clap(
        short = 'B',
        long = "buffer-size",
        about = "Size of reading buffer in bytes, default: 8192"
    )]
    pub buffer_size: Option<usize>,

    /// Number of thread crazyqc can use
    #[clap(
        short = 't',
        long = "threads",
        about = "Number of thread use by crazyqc, 0 use all avaible core, default: 0"
    )]
    pub threads: Option<usize>,

    /// Verbosity level
    #[clap(
        short = 'v',
        long = "verbosity",
        parse(from_occurrences),
        about = "verbosity level also control by environment variable CRAZYQC_LOG if flag is set CRAZYQC_LOG value is ignored"
    )]
    pub verbosity: i8,
}

/// Convert verbosity level (number of v) is log::Level
pub fn i82level(level: i8) -> Option<log::Level> {
    match level {
        std::i8::MIN..=0 => None,
        1 => Some(log::Level::Error),
        2 => Some(log::Level::Warn),
        3 => Some(log::Level::Info),
        4 => Some(log::Level::Debug),
        5..=std::i8::MAX => Some(log::Level::Trace),
    }
}

/// Set number of global rayon thread pool
pub fn set_nb_threads(nb_threads: usize) {
    log::info!("Set number of threads to {}", nb_threads);

    rayon::ThreadPoolBuilder::new()
        .num_threads(nb_threads)
        .build_global()
        .unwrap();
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn loglevel() {
        assert_eq!(i82level(i8::MIN), None);
        assert_eq!(i82level(-3), None);
        assert_eq!(i82level(1), Some(log::Level::Error));
        assert_eq!(i82level(2), Some(log::Level::Warn));
        assert_eq!(i82level(3), Some(log::Level::Info));
        assert_eq!(i82level(4), Some(log::Level::Debug));
        assert_eq!(i82level(5), Some(log::Level::Trace));
        assert_eq!(i82level(i8::MAX), Some(log::Level::Trace));
    }

    #[test]
    fn change_number_of_thread() {
        set_nb_threads(16);
        assert_eq!(rayon::current_num_threads(), 16);
    }
}
