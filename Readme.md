![Test](https://github.com/natir/crazyqc/workflows/Test/badge.svg)
![Lints](https://github.com/natir/crazyqc/workflows/Lints/badge.svg)
![MSRV](https://github.com/natir/crazyqc/workflows/MSRV/badge.svg)
[![CodeCov](https://codecov.io/gh/natir/crazyqc/branch/main/graph/badge.svg)](https://codecov.io/gh/natir/crazyqc)
[![Documentation](https://github.com/natir/crazyqc/workflows/Documentation/badge.svg)](https://natir.github.io/crazyqc/crazyqc)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/crazyqc/blob/main/LICENSE)

# CrazyQC 🧬 💻

CrazyQC is a try to generate some simple statistics on sequence data.

**WARNING**:
- CrazyQC has not yet been evaluated or even compared to *any* sequence QC tools
- CrazyQC is tested only on Linux
- CrazyQC is still in developpement many thing can change or be break

## Usage

```
crazyqc -q file1.fastq file2.fq ... zzz.fastq
```

By default CrazyQC use all avaible core you can control it with option `threads`:

```
crazyqc -t {number of thread} -q {your fastq file}
```

CrazyQC can also read bam file

```
crazyqc -t {number of thread} -q {your fastq file} -b {your bam file}
```


### Full usage

```
crazyqc 0.1
Pierre Marijon <pierre.marijon-ext@aphp.fr>
A crazily fast tools to compute small set of sequence quality control

USAGE:
    crazyqc [FLAGS] [OPTIONS]

FLAGS:
    -h, --help         Prints help information
    -v, --verbosity    verbosity level also control by environment variable CRAZYQC_LOG if flag is
                       set CRAZYQC_LOG value is ignored
    -V, --version      Prints version information

OPTIONS:
    -b, --bam <bam>...                 Bam input, optional
    -B, --buffer-size <buffer-size>    Size of reading buffer in bytes, default: 8192
    -q, --fastq <fastq>...             Fastq input
    -o, --output <output>              Path where result will be write, default: stdout
    -t, --threads <threads>            Number of thread use by crazyqc, 0 use all avaible core,
                                       default: 0
```

## Instalation

### With rust environment

If you haven't a rust environment you can use [rustup](https://rustup.rs/) or your package manager.

#### With cargo

```
cargo install --git https://github.com/natir/rustyread.git --tag 0.4
```

### From source

```
git clone https://github.com/natir/rustyread.git
cd rustyread
git checkout 0.4
cargo install --path .
```


## Minimum supported Rust version

Currently the minimum supported Rust version is 1.56.0.
