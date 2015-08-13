Duplex Novo
===========

This is a simple pipeline to process duplex sequencing data without the use of a reference sequence.

### Requirements

The pipeline requires a Unix command line and standard tools, and it must be able to find the `mafft` command on your `PATH`.

All known requirements are below. Version numbers in parentheses are what the development environment uses. Version numbers in **bold** are known to be required.

* [MAFFT](http://mafft.cbrc.jp/alignment/software/) (v7.123b)
* [Python](https://www.python.org/) (**2.7**)
* [awk](https://www.gnu.org/software/gawk/) (4.0.1)
* [gcc](https://gcc.gnu.org/) (4.8.4)

### Installation

Compile the C modules by `cd`ing to the top-level directory (where the file `Makefile` is) and run the command `make`.

### Usage
  
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
