//! Parse bam

/* std use */

/* crates use */
use anyhow::Context;

/* project use */
use crate::error;

/// Open a bam
fn open(
    path: &str,
    buffer_size: usize,
) -> anyhow::Result<noodles::bam::Reader<std::io::BufReader<std::fs::File>>> {
    log::debug!("Open file {}", path);

    let mut reader = noodles::bam::Reader::new(std::io::BufReader::with_capacity(
        buffer_size,
        std::fs::File::open(path)?,
    ));
    reader.read_header()?;
    reader.read_reference_sequences()?;

    Ok(reader)
}

/// Struct to parse Bam file
pub struct Bam {
    buffer_size: usize,
    paths: Vec<String>,
    local_record: noodles::bam::Record,
    current_input: noodles::bam::Reader<std::io::BufReader<std::fs::File>>,
}

impl Bam {
    /// Create a Bam struct, with inputs path and size of read buffer
    /// If last file of inputs can't be open this function return an anyhow::Error
    /// Other inputs can't be open is ignored, an error message is put in log
    pub fn new(mut inputs: Vec<String>, buffer_size: usize) -> anyhow::Result<Self> {
        let first_path = inputs
            .pop()
            .ok_or_else(|| anyhow::anyhow!("inputs is empty"))?;

        Ok(Self {
            buffer_size,
            paths: inputs,
            local_record: noodles::bam::Record::default(),
            current_input: open(&first_path, buffer_size).with_context(|| {
                error::Error::BamOpenError {
                    path: first_path.clone(),
                }
            })?,
        })
    }
}

impl Iterator for Bam {
    type Item = anyhow::Result<noodles::bam::record::Record>;

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
                            Err(e) => Some(Err(e.context(error::Error::BamOpenError {
                                path: new_path.clone(),
                            }))),
                        }
                    } else {
                        None // No new file end of iterator
                    }
                } else {
                    Some(Ok(self.local_record.clone()))
                }
            }
            Err(e) => Some(Err(
                anyhow::Error::new(e).context(error::Error::BamParsingError)
            )),
        }
    }
}

pub fn worker(wrapper: anyhow::Result<noodles::bam::record::Record>) -> (u64, u64, u64, u64) {
    match wrapper {
        Ok(seq) => {
            let mut at: u64 = 0;
            let mut gc: u64 = 0;
            let mut other: u64 = 0;

            for n in seq.sequence().bases() {
                match n {
                    noodles::bam::record::sequence::Base::A
                    | noodles::bam::record::sequence::Base::T => at += 1,
                    noodles::bam::record::sequence::Base::C
                    | noodles::bam::record::sequence::Base::G => gc += 1,
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

    fn create_bam_file() -> (NamedTempFile, String) {
        let tmp_file = tempfile::NamedTempFile::new().unwrap();

        {
            let writer = tmp_file.reopen().unwrap();
            let mut bam_writer = noodles::bam::Writer::new(writer);

            let header = noodles::sam::Header::builder().build();
            bam_writer.write_header(&header).unwrap();
            bam_writer
                .write_reference_sequences(header.reference_sequences())
                .unwrap();

            let record = noodles::sam::Record::builder()
                .set_read_name("1".parse().unwrap())
                .set_flags(noodles::sam::record::Flags::UNMAPPED)
                .set_sequence("ACTG".parse().unwrap())
                .build()
                .unwrap();
            bam_writer
                .write_sam_record(header.reference_sequences(), &record)
                .unwrap();

            let record = noodles::sam::Record::builder()
                .set_read_name("2".parse().unwrap())
                .set_flags(noodles::sam::record::Flags::UNMAPPED)
                .set_sequence("ACTGACTG".parse().unwrap())
                .build()
                .unwrap();
            bam_writer
                .write_sam_record(header.reference_sequences(), &record)
                .unwrap();

            let record = noodles::sam::Record::builder()
                .set_read_name("3".parse().unwrap())
                .set_flags(noodles::sam::record::Flags::UNMAPPED)
                .set_sequence("AACACGTGAGTCCGCACACCGGACG".parse().unwrap())
                .build()
                .unwrap();
            bam_writer
                .write_sam_record(header.reference_sequences(), &record)
                .unwrap();
        }

        let path = tmp_file.path().to_str().unwrap().to_string();

        (tmp_file, path)
    }

    #[test]
    fn open_bam_file() {
        let (tmp_file, path) = create_bam_file();

        assert!(open(&path, 10).is_ok());

        tmp_file.close().unwrap();

        assert!(open(&path, 10).is_err());
    }

    #[test]
    fn create_bam_parser() {
        let (tmp_file, path) = create_bam_file();

        assert!(Bam::new(vec![path.clone()], 10).is_ok());

        tmp_file.close().unwrap();

        assert!(Bam::new(vec![path], 10).is_err());
    }

    #[test]
    fn iterate_over_fasta() {
        let (_e, path) = create_bam_file();

        let mut reader = Bam::new(vec![path], 10).unwrap();

        assert_eq!(&reader.next().unwrap().unwrap().sequence()[..], [18, 132]);
        assert_eq!(
            &reader.next().unwrap().unwrap().sequence()[..],
            [18, 132, 18, 132]
        );
        assert_eq!(
            reader.next().unwrap().unwrap().sequence()[..],
            [17, 33, 36, 132, 20, 130, 36, 33, 33, 34, 68, 18, 64]
        );
    }

    #[test]
    fn run_worker() {
        let r1 = noodles::bam::Record::try_from_sam_record(
            noodles::sam::Header::builder()
                .build()
                .reference_sequences(),
            &noodles::sam::Record::builder()
                .set_read_name("1".parse().unwrap())
                .set_flags(noodles::sam::record::Flags::UNMAPPED)
                .set_sequence("ACTG".parse().unwrap())
                .build()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(worker(Ok(r1)), (2, 2, 0, 1));

        assert_eq!(worker(Err(anyhow::anyhow!("prout"))), (0, 0, 0, 0));
    }
}
