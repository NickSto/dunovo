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
AwkScript1='{print substr($2, 1, 12) substr($6, 1, 12), $1, substr($2, 18), $4, $5, substr($6, 18), $8}'
# This converts the barcode sequence into a canonical form that should be the same regardless of the
# order of the alpha and beta portions. It does this simply by sorting them based on their character
# values. So the output should be the same for a beta-alpha or an alpha-beta from the same fragment.
AwkScript2='BEGIN {
  FS = "\t"; OFS = "\t";
  BAR_LEN = 24;
  HALF = int(BAR_LEN/2);
  # Fill an array mapping each character to its ASCII value.
  for (i = 0; i < 128; i++) {
    char = sprintf("%c", i);
    ORD[char] = i;
  }
}
function get_sum(seq) {
  sum = 0;
  for (i = 1; i <= length(seq); i++) {
    sum += ORD[substr(seq, i, 1)];
  }
  return sum;
}
# Expect first and third column to be the FASTA header, like ">AAAAAAAACCAAAACCATCTCTTA.2".
length($1) == BAR_LEN + 3 && substr($1, 2, BAR_LEN) == substr($3, 2, BAR_LEN) {
  alpha = substr($1, 2, HALF);
  beta = substr($1, HALF+2, HALF);
  alpha_sum = get_sum(alpha);
  beta_sum = get_sum(beta);
  if (alpha_sum < beta_sum) {
    canonical = alpha beta;
    order = "ab";
  } else {
    canonical = beta alpha;
    order = "ba";
  }
  print canonical, order, $2, $4;
}'

function main {
  fastq1="$1"
  fastq2="$2"
  sscs="$3"
  # This transforms the input fastq's into a format that can be sorted by family with the "sort"
  # command. Mainly, it puts all the data for both read pairs on one line, and adds a column with
  # the barcode.
  # Warning: It assumes the fastq's have 4 lines per read!
  cat "$fastq1" | paste - - - - \
    | paste - <(cat "$fastq2" | paste - - - -) \
    | awk -F '\t' -v OFS='\t' "$AwkScript1" \
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
