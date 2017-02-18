#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Usage="Usage: \$ $(basename $0) [correct.py options] families.tsv > families.corrected.tsv"

function main {

  if [[ $# -lt 1 ]] || [[ $1 == '-h' ]]; then
    fail "$Usage"
  fi

  families="${@:$#:1}"

  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  script_dir=$(dirname "$script_path")

  bash "$script_dir/baralign.sh" "$families" refdir barcodes.bam
  samtools view -f 256 barcodes.bam \
    | python2 "$script_dir/correct.py" "$@" refdir/barcodes.fa \
    | sort
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
