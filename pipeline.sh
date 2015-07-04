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
# This translates the barcode sequence into a code that is independent of the order of the alpha and
# beta portions. So the output should be the same for a beta-alpha or an alpha-beta from the same
# fragment.
# Warning: This assumes all bases in the barcode are one of (uppercase) "ACGT". Characters not in
# that set will be ignored!
AwkScript2='BEGIN {
  BAR_LEN = 24;
  HALF = int(BAR_LEN/2);
  TR["AA"] = "A"; TR["AC"] = "B"; TR["AG"] = "D"; TR["AT"] = "E";
  TR["CA"] = "B"; TR["CC"] = "C"; TR["CG"] = "F"; TR["CT"] = "H";
  TR["GA"] = "D"; TR["GC"] = "F"; TR["GG"] = "G"; TR["GT"] = "I";
  TR["TA"] = "E"; TR["TC"] = "H"; TR["TG"] = "I"; TR["TT"] = "T";
}
# Expect first and third column to be the FASTA header, like ">AAAAAAAACCAAAACCATCTCTTA.2".
length($1) == BAR_LEN + 3 && substr($1, 2, BAR_LEN) == substr($3, 2, BAR_LEN) {
  code = ""
  for (i = 2; i <= HALF+1; i++) {
    j = i + HALF;
    pair = substr($1, i, 1) substr($1, j, 1);
    code = code TR[pair];
  }
  print code, substr($1, 2, BAR_LEN), $2, $4;
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
