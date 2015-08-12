#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

# At the moment this isn't really a production version of the pipeline.
# It mainly just documents how the commands are used.
function main {
  fastq1="$1"
  fastq2="$2"
  sscs="$3"
  # This transforms the input fastq's into a format that can be sorted by family with the "sort"
  # command. Mainly, it puts all the data for both read pairs on one line, and adds a column with
  # the barcode.
  # Warning: It assumes the fastq's have 4 lines per read!
  gunzip -c "$fastq1" | paste - - - - \
    | paste - <(gunzip -c "$fastq2" | paste - - - -) \
    | awk -f make-barcodes.awk \
    | sort \
    | align_families.py \
    | duplex.py \
    | gzip -c - \
    > "$sscs"
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
