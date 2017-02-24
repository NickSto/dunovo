#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

MinReadsDefault=3
QualThresDefault=25
QualFormatDefault=sanger
Usage="Usage: \$ $(basename $0) [options] families.msa.tsv dcs1.fa dcs2.fa [sscs1.fa sscs2.fa]"

function main {

  # Read arguments.
  if [[ $# -lt 2 ]]; then
    fail "$Usage"
  fi
  min_reads=$MinReadsDefault
  qual_thres=$QualThresDefault
  qual_format=$QualFormatDefault
  while getopts ":r:q:F:h" opt; do
  case "$opt" in
      r) min_reads=$OPTARG;;
      q) qual_thres=$OPTARG;;
      F) qual_format=$OPTARG;;
      h) echo "$USAGE"
         exit;;
    esac
  done
  # Get positionals.
  alignments="${@:$OPTIND:1}"
  dcs1="${@:$OPTIND+1:1}"
  dcs2="${@:$OPTIND+2:1}"
  sscs1="${@:$OPTIND+3:1}"
  sscs2="${@:$OPTIND+4:1}"

  keep_sscs=
  if [[ $sscs1 ]] && [[ $sscs2 ]]; then
    keep_sscs=true
  elif [[ $sscs1 ]] || [[ $sscs2 ]]; then
    fail "Error: Must give both sscs1 and sscs2 or neither."
  fi

  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  script_dir=$(dirname "$script_path")

  sscs_args=
  if [[ $keep_sscs ]]; then
    sscs_args='--sscs-file sscs.fa'
  fi

  # Locate outconv.py.
  # In $script_dir/utils in a normal installation, $script_dir in a Conda installation.
  outconv_script="$script_dir/utils/outconv.py"
  if ! [[ -f "$outconv_script" ]]; then
    outconv_script="$script_dir/outconv.py"
  fi

  python2 "$script_dir/dunovo.py" -r $min_reads -q $qual_thres -F $qual_format "$alignments" \
    $sscs_args > duplex.fa
  python2 "$outconv_script" duplex.fa -1 "$dcs1" -2 "$dcs2"

  if [[ $keep_sscs ]]; then
    python2 "$outconv_script" sscs.fa -1 "$sscs1" -2 "$sscs2"
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
