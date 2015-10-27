#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

BarcodeLen=12
SpacerLen=5
ReadlenDefault=101
BwaCmd=${BwaCmd:="bwa"}

Usage="Usage: \$ $(basename $0) ref.fa reads_1.fq reads_2.fq [readlen]
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
  elif [[ $# -lt 3 ]] || [[ $1 == '-h' ]]; then
    fail "$Usage"
  else
    ref="$1"
    fastq1="$2"
    fastq2="$3"
  fi
  if [[ $# -ge 4 ]]; then
    readlen="$4"
  fi

  print_versions $script_dir

  rlenreal=$((readlen-BarcodeLen-SpacerLen))

  python $script_dir/tag_to_header.py --infile1 $fastq1 --infile2 $fastq2 --outfile1 read_1.fq.smi \
    --outfile2 read_2.fq.smi --barcode_length $BarcodeLen --spacer_length $SpacerLen
  $BwaCmd aln $ref read_1.fq.smi > read_1.aln 
  $BwaCmd aln $ref read_2.fq.smi > read_2.aln
  $BwaCmd sampe -s $ref read_1.aln read_2.aln read_1.fq.smi read_2.fq.smi > PE_reads.sam
  samtools view -Sbu PE_reads.sam | samtools sort - PE_reads.sort
  python $script_dir/ConsensusMaker.py --infile PE_reads.sort.bam --tagfile PE_reads.tagcounts \
    --outfile SSCS.bam --min 3 --max 1000 --cutoff 0.7 --Ncutoff 0.3 --readlength $rlenreal \
    --read_type dpm --filt osn
  samtools view -bu SSCS.bam | samtools sort - SSCS.sort
  python $script_dir/DuplexMaker.py --infile SSCS.sort.bam --outfile DCS_data --Ncutoff 0.3 \
    --readlength $rlenreal
  $BwaCmd aln $ref DCS_data.r1.fq > DCS_data.r1.aln
  $BwaCmd aln $ref DCS_data.r2.fq > DCS_data.r2.aln
  $BwaCmd sampe -s $ref DCS_data.r1.aln DCS_data.r2.aln DCS_data.r1.fq DCS_data.r2.fq \
    > DCS_PE.aln.sam
  samtools view -Sbu DCS_PE.aln.sam | samtools sort - DCS_PE.aln.sort
  samtools index DCS_PE.aln.sort.bam
  samtools view -F 4 -b DCS_PE.aln.sort.bam > DCS_PE.filt.bam
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
  python --version 2>&1 | sed -E 's/python\s//I'
  echo -en 'BWA:\t\t0.6.2\t\t'
  $BwaCmd 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  echo -en 'Samtools:\t0.1.18\t\t'
  samtools 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  echo -en 'PySAM:\t\t0.7.5\t\t'
  python -c 'import pysam; print pysam.__version__'
  echo -en 'BioPython:\t1.62\t\t'
  python -c 'import Bio; print Bio.__version__'
}


function fail {
  echo "$@" >&2
  exit 1
}


main "$@"
