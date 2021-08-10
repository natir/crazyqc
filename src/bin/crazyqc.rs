/* crate use */
use clap::Clap;
use rayon::prelude::*;

/* local use */
use crazyqc::*;

#[cfg(not(tarpaulin_include))]
fn main() -> anyhow::Result<()> {
    let params = cli::Command::parse();

    /* Manage cli */
    if let Some(level) = cli::i82level(params.verbosity) {
        env_logger::builder()
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .filter_level(level.to_level_filter())
            .init();
    } else {
        env_logger::Builder::from_env("CRAZYQC_LOG")
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .init();
    }

    if let Some(threads) = params.threads {
        cli::set_nb_threads(threads);
    }

    let buffer_size: usize = if let Some(b) = params.buffer_size {
        b
    } else {
        8192
    };

    let mut output: Box<dyn std::io::Write> = if let Some(o) = params.output {
        Box::new(std::io::BufWriter::new(std::fs::File::open(o)?))
    } else {
        Box::new(std::io::stdout())
    };

    log::info!("Start read fastq");
    let reader = input::Fastq::new(params.fastq, buffer_size)?;

    let (at, gc, other, nb_reads) = reader
        .par_bridge()
        .map(input::fastq::worker)
        .reduce(reduce_identity, reduce);

    let total = (at + gc + other) as f64;
    log::info!("End read fastq");

    writeln!(output, "type,at,gc,other,mean_length,n")?;
    writeln!(
        output,
        "fastq,{:.4},{:.4},{:.4},{:.2},{}",
        (at as f64 / total) * 100.0,
        (gc as f64 / total) * 100.0,
        (other as f64 / total) * 100.0,
        total / nb_reads as f64,
        nb_reads
    )?;

    /* Run count of bam file if option is set */
    if let Some(bams_path) = params.bam {
        log::info!("Start read bam");
        let reader = input::Bam::new(bams_path, buffer_size)?;

        let (at, gc, other, nb_reads) = reader
            .par_bridge()
            .map(input::bam::worker)
            .reduce(reduce_identity, reduce);
        log::info!("End read bam");

        let total = (at + gc + other) as f64;

        writeln!(
            output,
            "bam,{:.4},{:.4},{:.4},{:.2},{}",
            at as f64 / total * 100.0,
            gc as f64 / total * 100.0,
            other as f64 / total * 100.0,
            total / nb_reads as f64,
            nb_reads
        )?;
    }

    Ok(())
}

fn reduce(a: (u64, u64, u64, u64), b: (u64, u64, u64, u64)) -> (u64, u64, u64, u64) {
    (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3)
}

fn reduce_identity() -> (u64, u64, u64, u64) {
    (0, 0, 0, 0)
}

#[cfg(test)]
mod t {
    use super::*;

    #[test]
    fn reduce_() {
        assert_eq!(reduce((1, 1, 1, 1), (2, 2, 2, 2)), (3, 3, 3, 3));
    }

    #[test]
    fn reduce_identity_() {
        assert_eq!(reduce_identity(), (0, 0, 0, 0))
    }
}
