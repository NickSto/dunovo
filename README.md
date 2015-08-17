Duplex Novo
===========

This is a simple pipeline to process duplex sequencing data without the use of a reference sequence.


### Requirements

The pipeline requires a Unix command line, and it must be able to find the `mafft` command on your [`PATH`](https://en.wikipedia.org/wiki/Search_path).

All known requirements are below. Version numbers in parentheses are what the development environment uses. Version numbers in **bold** are known to be required.

* [MAFFT](http://mafft.cbrc.jp/alignment/software/) (v7.123b)
* [Python](https://www.python.org/) (**2.7**)
* [gcc](https://gcc.gnu.org/) (4.8.4)
Standard unix tools:
* [bash](https://www.gnu.org/software/bash/bash.html) (4.0)
* [awk](https://www.gnu.org/software/gawk/) (4.0.1)
* [paste](https://www.gnu.org/software/coreutils/coreutils.html), [sort](https://www.gnu.org/software/coreutils/coreutils.html), [cat](https://www.gnu.org/software/coreutils/coreutils.html) (8.21)


### Installation

`git clone` the source to any directory, or click "Download ZIP", unzip it, and place the "duplex-master" directory anywhere.

You'll need to compile the C modules before using it. Do this in a terminal by `cd`ing to the source directory (where the file `Makefile` is) and run the command `make`.


### Usage

This example shows how to go from raw duplex sequencing data to the final duplex consensus sequences.

The example assumes you want to process duplex reads in the files `reads_1.fastq` and `reads_2.fastq`, and have `cd`'d to their directory. It also assumes you've placed the commands `align_families.py` and `duplex.py` on your `PATH`. Note that where it says `make-barcodes.awk`, you should replace that with the actual path to the script `make-barcodes.awk` included in this pipeline.

1. Sort the reads into families based on their barcodes.  
    ```bash
    $ cat reads_1.fastq | paste - - - - \
      | paste - <(cat reads_2.fastq | paste - - - -) \
      | awk -f make-barcodes.awk \
      | sort > families.tsv
    ```

2. Do multiple sequence alignments of the read families.  
`$ align_families.py families.tsv > families.msa.tsv`

3. Build duplex consensus sequences of from the aligned families.  
`$ duplex.py families.msa.tsv > duplex.fa`

See all options for a given command by giving it the `-h` flag.


### Details

##### 1. Sort the reads into families based on their barcodes.  

    $ cat reads_1.fastq | paste - - - - \
      | paste - <(cat reads_2.fastq | paste - - - -) \
      | awk -f make-barcodes.awk \
      | sort > families.tsv

This command pipeline will transform each pair of reads into a one-line record, split the 12bp barcodes off them, and sort by their combined barcode. The end result is a file (named `families.tsv` above) listing read pairs, grouped by barcode. See `make-barcodes.awk` for the details on the formation of the barcodes and the format.

Note: This step requires your FASTQ files to have exactly 4 lines per read (no multi-line sequences). Also, in the output, the read sequence does not include the barcode or the 5bp constant sequence after it. You can customize the length of the barcode or constant sequence by setting the awk constants `BAR_LEN` and `INVARIANT` (i.e. `awk -v BAR_LEN=10 make-barcodes.awk`).


##### 2. Do multiple sequence alignments of the read families.  

`$ align_families.py families.tsv > families.msa.tsv`

This step aligns each family of reads, but it processes each strand separately. It can be parallelized with the `-p` option, but at the moment that will cause the output to only be generated at the end, instead of streaming it as it's generated.


##### 3. Build duplex consensus sequences of from the aligned families.  

`$ duplex.py families.msa.tsv > duplex.fa`

This calls a consensus sequence from the multiple sequence alignments of the previous step. It does this in two steps: First, single-strand consensus sequences (SSCSs) are called from the family alignments, then duplex consensus sequences are called from pairs of SSCSs.

When calling SSCSs, by default 3 reads are required to successfully create a consensus from each strand. Quality filtering is done at this step by excluding bases below a quality threshold. By default, no base with a PHRED quality less than 20 will contribute to the consensus. If no base passes the threshold or there is no majority base, an `N` will be inserted.

The duplex consensus sequences are created by comparing the two SSCSs. For each base, if they agree, that base will be inserted. If they disagree, the IUPAC ambiguity code for the two bases will be used. Note that a disagreement between a base and a gap will result in an `N`. A planned feature is to use information from the raw reads contributing to the duplex to make a call in such a case, coding uncertainty into quality scores.

The output of this step is the duplex consensus sequences in FASTA format. By default, it will only include full duplex consensuses, meaning if one of the two SSCSs are missing, that sequence will be omitted. But these sequences can be included with the `--incl-sscs` option, which will add lone SSCSs to the output.
