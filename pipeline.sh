#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

# The awk code that transforms the one-line fastq record pair into the output that can be sorted
# by barcode. Column 1 is the barcode, which is simply the first 12bp of read 1 concatenated with
# the first 12bp of read 2. Columns 3 and 6 are the sequences of reads 1 and 2; the 12bp barcode and
# 5bp constant sequence are removed.
AwkScript='{print substr($2, 1, 12) substr($6, 1, 12), $1, substr($2, 18), $4, $5, substr($6, 18), $8}'

function main {
  fastq1="$1"
  fastq2="$2"
  sscs="$3"
  cat "$fastq1" | paste - - - - \
    | paste - <(cat "$fastq2" | paste - - - -) \
    | awk -F '\t' -v OFS='\t' "$AwkScript" \
    | sort -k 1 \
    | sscs.py \
    | gzip -c - \
    > "$sscs"
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
