//! Parse fastq

/* std use */

/* crates use */
use anyhow::Context;

/* project use */
use crate::error;

/// Open a fasta
fn open(
    path: &str,
    buffer_size: usize,
) -> anyhow::Result<noodles::fastq::Reader<std::io::BufReader<Box<dyn std::io::Read + Send>>>> {
    log::debug!("Open file {}", path);

    Ok(noodles::fastq::Reader::new(
        std::io::BufReader::with_capacity(
            buffer_size,
            niffler::send::get_reader(Box::new(std::fs::File::open(path)?))?.0,
        ),
    ))
}

/// Struct to parse Fastq file
pub struct Fastq {
    buffer_size: usize,
    paths: Vec<String>,
    local_record: noodles::fastq::Record,
    current_input: noodles::fastq::Reader<std::io::BufReader<Box<dyn std::io::Read + Send>>>,
}

impl Fastq {
    /// Create a Fastq struct, with inputs path and size of read buffer
    /// If last file of inputs can't be open this function return an anyhow::Error
    /// Other inputs can't be open is ignored, an error message is put in log
    pub fn new(mut inputs: Vec<String>, buffer_size: usize) -> anyhow::Result<Self> {
        let first_path = inputs
            .pop()
            .ok_or_else(|| anyhow::anyhow!("inputs is empty"))?;

        Ok(Self {
            buffer_size,
            paths: inputs,
            local_record: noodles::fastq::Record::default(),
            current_input: open(&first_path, buffer_size).with_context(|| {
                error::Error::FastqOpenError {
                    path: first_path.clone(),
                }
            })?,
        })
    }
}

impl Iterator for Fastq {
    type Item = anyhow::Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.current_input.read_record(&mut self.local_record) {
            Ok(nb_bytes_read) => {
                if nb_bytes_read == 0 {
                    if let Some(new_path) = self.paths.pop() {
                        match open(&new_path, self.buffer_size) {
                            Ok(new_input) => {
                                self.current_input = new_input;
                                self.next()
                            }
                            Err(e) => Some(Err(e.context(error::Error::FastqOpenError {
                                path: new_path.clone(),
                            }))),
                        }
                    } else {
                        None // No new file end of iterator
                    }
                } else {
                    Some(Ok(self.local_record.sequence().to_vec()))
                }
            }
            Err(e) => Some(Err(
                anyhow::Error::new(e).context(error::Error::FastqParsingError)
            )),
        }
    }
}

pub fn worker(wrapper: anyhow::Result<Vec<u8>>) -> (u64, u64, u64, u64) {
    match wrapper {
        Ok(seq) => {
            let mut at: u64 = 0;
            let mut gc: u64 = 0;
            let mut other: u64 = 0;

            for n in seq {
                match n {
                    b'a' | b'A' | b't' | b'T' => at += 1,
                    b'c' | b'C' | b'g' | b'G' => gc += 1,
                    _ => other += 1,
                }
            }

            (at, gc, other, 1)
        }
        Err(e) => {
            e.chain().for_each(|b| log::error!("{}", b));
            (0, 0, 0, 0)
        }
    }
}

#[cfg(test)]
mod t {
    use tempfile::NamedTempFile;

    use super::*;
    use std::io::Write;

    fn create_fastq_file() -> (NamedTempFile, String) {
        let tmp_file = tempfile::NamedTempFile::new().unwrap();

        {
            let mut writer = tmp_file.reopen().unwrap();
            writeln!(writer, "@1\nACTG\n+\n!!!!").unwrap();
            writeln!(writer, "@2\nACTGACTG\n+\n!!!!!!!!").unwrap();
            write!(
                writer,
                "@3\nAACACGTGAGTCCGCACACCGGACG\n+\n`kL7sm$xKvE8.RT`[kgO!34O'"
            )
            .unwrap();
        }

        let path = tmp_file.path().to_str().unwrap().to_string();

        (tmp_file, path)
    }

    #[test]
    fn open_fastq_file() {
        let (tmp_file, path) = create_fastq_file();

        assert!(open(&path, 10).is_ok());

        tmp_file.close().unwrap();

        assert!(open(&path, 10).is_err());
    }

    #[test]
    fn create_fastq_parser() {
        let (tmp_file, path) = create_fastq_file();

        assert!(Fastq::new(vec![path.clone()], 10).is_ok());

        tmp_file.close().unwrap();

        assert!(Fastq::new(vec![path], 10).is_err());
    }

    #[test]
    fn iterate_over_fasta() {
        let (_e, path) = create_fastq_file();

        let mut reader = Fastq::new(vec![path], 10).unwrap();
        assert_eq!(reader.next().unwrap().unwrap(), b"ACTG".to_vec());
        assert_eq!(reader.next().unwrap().unwrap(), b"ACTGACTG".to_vec());
        assert_eq!(
            reader.next().unwrap().unwrap(),
            b"AACACGTGAGTCCGCACACCGGACG".to_vec()
        );
    }

    #[test]
    fn run_worker() {
        assert_eq!(worker(Ok(b"ACTG".to_vec())), (2, 2, 0, 1));

        assert_eq!(worker(Err(anyhow::anyhow!("prout"))), (0, 0, 0, 0));
    }
}
