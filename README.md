# _Du Novo_

This is a pipeline for processing of duplex sequencing data without the use of a reference genome.

The pipeline was designed for use with the duplex method described in [Kennedy *et al.* 2014](https://dx.doi.org/10.1038/nprot.2014.170), but the assumptions are relatively minimal, so you should be able to apply it to variants of the protocol.

## Using _Du Novo_

_Du Novo_ can be used in one of two ways:

 * via Galaxy
 * on the command line

## Running _Du Novo_ from Galaxy

We created a comprehensive [tutorial](https://github.com/galaxyproject/dunovo/wiki) explaining all aspects of interactive use of _De Novo_ from within [Galaxy](http://usegalaxy.org).

## Running _Du Novo_ on the command line

### Requirements

The pipeline requires a Unix command line, and it must be able to find the `mafft` command on your [`PATH`](https://en.wikipedia.org/wiki/Search_path).

All known requirements are below. Version numbers in parentheses are what the development environment uses. Version numbers in **bold** are known to be required.

* [MAFFT](http://mafft.cbrc.jp/alignment/software/) (v7.123b)
* [Python](https://www.python.org/) (**2.7**)  
* And standard unix tools:
 -  [gcc](https://gcc.gnu.org/) (4.8.4)
 -  [make](https://www.gnu.org/software/make/) (3.81)
 -  [bash](https://www.gnu.org/software/bash/bash.html) (4.0)
 -  [awk](https://www.gnu.org/software/gawk/) (4.0.1)
 -  [paste](https://www.gnu.org/software/coreutils/coreutils.html) (8.21)
 -  [sort](https://www.gnu.org/software/coreutils/coreutils.html) (8.21)

### Installation

    $ git clone git@github.com:makrutenko/duplex.git
    $ cd duplex
    $ make

Instead of the `git` command, you can just click "Download ZIP", unzip it, and `cd` to the "duplex-master" directory.

You'll need to compile the C modules before using it. Do this in a terminal by `cd`ing to the source directory (where the file `Makefile` is) and run the command `make`.


### Usage

This example shows how to go from raw duplex sequencing data to the final duplex consensus sequences.

Your raw reads should be in `reads_1.fastq` and `reads_2.fastq`. And the scripts `align_families.py` and `dunovo.py` should be on your `PATH`. Also, in the following command, replace `make-barcodes.awk` with the actual path to that script (included in this pipeline).

1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  
    ```bash
    $ paste reads_1.fastq reads_2.fastq \
      | paste - - - - \
      | awk -f make-barcodes.awk \
      | sort > families.tsv
    ```

2. Do multiple sequence alignments of the read families.  
`$ align_families.py families.tsv > families.msa.tsv`

3. Build duplex consensus sequences from the aligned families.  
`$ dunovo.py families.msa.tsv > duplex.fa`

See all options for a given command by giving it the `-h` flag.


### Details

#### 1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  

    $ paste reads_1.fastq reads_2.fastq \
      | paste - - - - \
      | awk -f make-barcodes.awk \
      | sort > families.tsv

This command pipeline will transform each pair of reads into a one-line record, split the 12bp barcodes off them, and sort by their combined barcode. The end result is a file (named `families.tsv` above) listing read pairs, grouped by barcode. See `make-barcodes.awk` for the details on the formation of the barcodes and the format.

Note: This step requires your FASTQ files to have exactly 4 lines per read (no multi-line sequences). Also, in the output, the read sequence does not include the barcode or the 5bp constant sequence after it. You can customize the length of the barcode or constant sequence by setting the awk constants `TAG_LEN` and `INVARIANT` (i.e. `awk -v TAG_LEN=10 make-barcodes.awk`).


#### 2. Do multiple sequence alignments of the read families.  

`$ align_families.py families.tsv > families.msa.tsv`

This step aligns each family of reads, but it processes each strand separately. It can be parallelized with the `-p` option.


#### 3. Build duplex consensus sequences from the aligned families.  

`$ dunovo.py families.msa.tsv > duplex.fa`

This calls a consensus sequence from the multiple sequence alignments of the previous step. It does this in two steps: First, single-strand consensus sequences (SSCSs) are called from the family alignments, then duplex consensus sequences are called from pairs of SSCSs.

When calling SSCSs, by default 3 reads are required to successfully create a consensus from each strand (change this with `-r`). Quality filtering is done at this step by excluding bases below a quality threshold. By default, no base with a PHRED quality less than 20 will contribute to the consensus (change this with `-q`). If no base passes the threshold or there is no majority base, `N` will be used.

The duplex consensus sequences are created by comparing the two SSCSs. For each base, if they agree, that base will be inserted. If they disagree, the IUPAC ambiguity code for the two bases will be used. Note that a disagreement between a base and a gap will result in an `N`.

The output of this step is the duplex consensus sequences in FASTA format. By default, it will only include full duplex consensuses, meaning if one of the two SSCSs are missing, that sequence will be omitted. Include these with the `--incl-sscs` option.

The reads will be printed in one, interleaved file, with the naming format:  
`>{barcode}.{mate} {# reads in strand 1 family}/{# reads in strand 2 family}`  
e.g.  
`>TTGCGCCAGGGCGAGGAAAATACT.1 8/13`

But this isn't easy to work with. A better output is in development, but for now you can use the awk script `outconv.awk` to convert the interleaved output file into two standard forward/reverse paired files with a standard naming convention:

    $ awk -f utils/outconv.awk -v target=1 duplex.fa > duplex_1.fa
    $ awk -f utils/outconv.awk -v target=2 duplex.fa > duplex_2.fa
