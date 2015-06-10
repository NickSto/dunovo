#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

AwkScript='{print substr($2, 1, 12) substr($6, 1, 12), $1, substr($2, 13), $4, $5, substr($6, 13), $8}'

function main {
  fastq1="$1"
  fastq2="$2"
  cat "$fastq1" | paste - - - - \
    | paste - <(cat "$fastq2" | paste - - - -) \
    | awk -F '\t' -v OFS='\t' "$AwkScript" \
    | sort -k 1 \
    | sscs.py
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
