//! All stuff relate to command line

/// Struct use to parse command line argument
#[derive(clap::Parser, Debug)]
#[clap(
    name = "crazyqc",
    version = "0.1 ",
    author = "Pierre Marijon <pierre.marijon-ext@aphp.fr>",
    override_help = "A crazily fast tools to compute small set of sequence quality control"
)]
pub struct Command {
    /// Fastq input
    #[clap(short = 'q', long = "fastq")]
    pub fastq: Vec<String>,

    /// Bam input
    #[clap(short = 'b', long = "bam")]
    pub bam: Option<Vec<String>>,

    /// Path where result will be write, default: stdout
    #[clap(short = 'o', long = "output")]
    pub output: Option<String>,

    /// Size of reading buffer in bytes, default: 8192
    #[clap(short = 'B', long = "buffer-size")]
    pub buffer_size: Option<usize>,

    /// Number of thread use by crazyqc, 0 use all avaible core, default: 0
    #[clap(short = 't', long = "threads")]
    pub threads: Option<usize>,

    /// Verbosity level also control by environment variable CRAZYQC_LOG if flag is set CRAZYQC_LOG value is ignored
    #[clap(short = 'v', long = "verbosity", parse(from_occurrences))]
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
