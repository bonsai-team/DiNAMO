
# DiNAMO: An exact and efficient method for IUPAC motif discovery in DNA sequences

[![Build Status](https://travis-ci.org/bonsai-team/DiNAMO.svg?branch=master)](https://travis-ci.org/bonsai-team/DiNAMO) [![Build status](https://ci.appveyor.com/api/projects/status/a9powubl5mqeigvt/branch/master?svg=true)](https://ci.appveyor.com/project/chadisaad/dinamo/branch/master)

The DiNAMO software implements an exhaustive algorithm to detect over-represented IUPAC motifs in a set of DNA sequences. It has two modes: *scanning* mode (default), where all windows are parsed, or *fixed-position* mode (see optional parameter *-p* ), where only motifs occurring at a specific position in the sequences are taken into account.  

DiNAMO can be used in a variety of applications, such as ChIP-seq peak analysis for transcription factor binding sites identification, or finding motifs that induce systematic sequencing errors.

## Installation

### Binaries
Binaries for Windows, OS X and Linux are available [here](https://github.com/bonsai-team/DiNAMO/releases).

### Build from source code
#### Dependencies

    Boost library (if not in the standard includes, please use -I /path/to/boost in the Makefile)


#### Compilation
    make

The executable file is located in the bin/ directory

## Usage

DiNAMO takes as input two fasta files, which contain respectively the positive dataset (signal.fa) and the negative dataset (control.fa). You should specify also the motif length '-l'. Optionally, you can provide a maximum number of degenerate letters '-d'. The output file contains the over-represented motifs found in the positive dataset compared to the negative dataset, in the MEME format (see http://meme-suite.org/doc/meme-format.html).

        dinamo -pf signal.fa  -nf control.fa -l 6

### Optional Parameters


* maximum number *i* of degenerate letters

        -d <i>

*  Only process motifs at a offset i (0..N)  related to the **end** of each sequence

        -p <i>

* output file

        -o <file>

* change Fisher's exact test p-value threshold *f* (default: 0.05)

        -t <f>

* do not display logs

        --no-log


## References

how to cite this tool:

>C. Saad, L. Noé, H. Richard, J. Leclerc, M.-P. Buisine, H. Touzet, and M. Figeac. Dinamo: highly sensitive dna motif discovery in high-throughput sequencing data. BMC Bioinformatics, 19(1):223, Jun 2018. [https://doi.org/10.1186/s12859-018-2215-1](https://doi.org/10.1186/s12859-018-2215-1)