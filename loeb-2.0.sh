#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

BarcodeLen=${BarcodeLen:=12}
SpacerLen=${SpacerLen:=5}
BwaCmd=${BwaCmd:="bwa"}
SamtoolsCmd=${SamtoolsCmd:="samtools"}
PythonCmd=${PythonCmd:="python"}
ReadlenDefault=101

Usage="Usage: \$ $(basename $0) ref.fa reads_1.fq reads_2.fq readlen [outdir]
Run the Loeb pipeline as it was published in the Kennedy et al. 2014 paper (release 2.0).
If readlen is not provided, it will assume ${ReadlenDefault}bp.
Dependencies:
Python >= 2.7 (and < 3.0)
BWA <= 0.6.2
Samtools <= 0.1.18
BioPython 1.62
PySAM 0.7.5
To just check your dependency versions, run \$ $(basename $0) -v"

function main {

  script_dir=$(real_dir)

  if [[ $# == 1 ]] && [[ $1 == '-v' ]]; then
    print_versions $script_dir
    exit
  elif [[ $# -lt 4 ]] || [[ $1 == '-h' ]]; then
    fail "$Usage"
  else
    ref="$1"
    fastq1="$2"
    fastq2="$3"
    readlen="$4"
  fi
  if [[ $# -ge 5 ]]; then
    outdir="$5"
  else
    outdir=.
  fi

  print_versions $script_dir

  echo "
Parameters:
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
readlen: $readlen
"

  rlenreal=$((readlen-BarcodeLen-SpacerLen))

  $PythonCmd $script_dir/tag_to_header.py --infile1 $fastq1 --infile2 $fastq2 \
    --outfile1 $outdir/read_1.fq.smi --outfile2 $outdir/read_2.fq.smi --barcode_length $BarcodeLen \
    --spacer_length $SpacerLen
  $BwaCmd aln $ref $outdir/read_1.fq.smi > $outdir/read_1.aln
  $BwaCmd aln $ref $outdir/read_2.fq.smi > $outdir/read_2.aln
  $BwaCmd sampe -s $ref $outdir/read_1.aln $outdir/read_2.aln \
    $outdir/read_1.fq.smi $outdir/read_2.fq.smi > $outdir/PE_reads.sam
  $SamtoolsCmd view -Sbu $outdir/PE_reads.sam | $SamtoolsCmd sort - $outdir/PE_reads.sort
  $PythonCmd $script_dir/ConsensusMaker.py --tagfile $outdir/PE_reads.tagcounts \
    --infile $outdir/PE_reads.sort.bam --outfile $outdir/SSCS.bam --min 3 --max 1000 --cutoff 0.7 \
    --Ncutoff 0.3 --readlength $rlenreal --read_type dpm --filt osn
  $SamtoolsCmd view -bu $outdir/SSCS.bam | $SamtoolsCmd sort - $outdir/SSCS.sort
  $PythonCmd $script_dir/DuplexMaker.py --infile $outdir/SSCS.sort.bam --outfile $outdir/DCS_data \
    --Ncutoff 0.3 --readlength $rlenreal
  $BwaCmd aln $ref $outdir/DCS_data.r1.fq > $outdir/DCS_data.r1.aln
  $BwaCmd aln $ref $outdir/DCS_data.r2.fq > $outdir/DCS_data.r2.aln
  $BwaCmd sampe -s $ref $outdir/DCS_data.r1.aln $outdir/DCS_data.r2.aln \
    $outdir/DCS_data.r1.fq $outdir/DCS_data.r2.fq > $outdir/DCS_PE.aln.sam
  $SamtoolsCmd view -Sbu $outdir/DCS_PE.aln.sam | $SamtoolsCmd sort - $outdir/DCS_PE.aln.sort
  $SamtoolsCmd index $outdir/DCS_PE.aln.sort.bam
  $SamtoolsCmd view -F 4 -b $outdir/DCS_PE.aln.sort.bam > $outdir/DCS_PE.filt.bam
}


# Get the script's actual directory path
function real_dir {
  # Does readlink -f work? (It doesn't on BSD.)
  if readlink -f dummy >/dev/null 2>/dev/null; then
    dirname $(readlink -f ${BASH_SOURCE[0]})
  else
    # If readlink -f doesn't work (like on BSD).
    # Read the link destination from the output of ls -l and cd to it.
    # Have to cd to the link's directory first, to handle relative links.
    # With help from https://stackoverflow.com/a/246128/726773
    unset CDPATH
    local source="${BASH_SOURCE[0]}"
    while [[ -h "$source" ]]; do
      local dir="$(cd -P $(dirname "$source") && pwd)"
      local link="$(ls -l "$source" | awk '{print $NF}')"
      # absolute or relative path?
      if [[ "$link" == /* ]]; then
        source="$link"
      else
        source="$dir/$link"
      fi
    done
    dir="$(cd -P $(dirname "$source") && pwd)"
    echo "$dir"
  fi
}


function print_versions {
  script_dir="$1"
  echo -e 'VERSIONS\trecommended\tpresent'
  echo -en 'pipeline:\te0897da\t\t'
  unset CDPATH
  cd $script_dir
  git log --oneline -n 1 | grep --color=never -Eo '^\S+'
  cd - >/dev/null
  echo -en 'Python:\t\t2.7\t\t'
  $PythonCmd --version 2>&1 | sed -E 's/python\s//I'
  echo -en 'BWA:\t\t0.6.2\t\t'
  $BwaCmd 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  echo -en 'Samtools:\t0.1.18\t\t'
  $SamtoolsCmd 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  echo -en 'PySAM:\t\t0.7.5\t\t'
  $PythonCmd -c 'import pysam; print pysam.__version__'
  echo -en 'BioPython:\t1.62\t\t'
  $PythonCmd -c 'import Bio; print Bio.__version__'
}


function fail {
  echo "$@" >&2
  exit 1
}


main "$@"
