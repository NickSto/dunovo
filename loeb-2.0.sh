#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Start=${Start:=}
BarcodeLen=${BarcodeLen:=12}
SpacerLen=${SpacerLen:=5}
StartClip=${StartClip:=5}
EndClip=${EndClip:=5}
BwaCmd=${BwaCmd:="bwa"}
SamtoolsCmd=${SamtoolsCmd:="samtools"}
PythonCmd=${PythonCmd:="python"}
JavaCmd=${JavaCmd:="java"}
PicardDir=${PicardDir:="$HOME/src/picard-tools-1.100"}
GatkDir=${GatkDir:="$HOME/src/GenomeAnalysisTK"}

Usage="Usage: \$ $(basename $0) [-d|-c|-a] ref.fa reads_1.fq reads_2.fq readlen [outdir]
Run the Loeb pipeline as it was published in the Kennedy et al. 2014 paper
(release 2.0). If -d (\"duplex\") is given, it will stop after producing the
final duplex reads (step 62). This is the default. If -c (\"cleanup\") is given,
it will skip producing the duplex reads, assuming it's already been done, and
just do the filtering, realignment, and trimming (steps 63-71). If -a (\"all\")
is given, it will do the whole pipeline (both halves). If it's not doing the
second part, Picard and GATK are not required. Otherwise, provide the paths to
the directories containing their .jars by setting \$PicardDir and \$GatkDir.
Dependencies:
Python    >= 2.7 (and < 3.0)
BWA       <= 0.6.2
Samtools  <= 0.1.18
BioPython 1.62
PySAM     0.7.5
Picard    1.107
GATK      2.4-9
To just check your dependency versions, run \$ $(basename $0) -v"

function main {

  script_dir=$(real_dir)

  duplex=true
  cleanup=''
  if [[ $# -ge 1 ]]; then
    if [[ $1 == '-v' ]]; then
      print_versions $script_dir
      exit
    elif [[ $1 == '-d' ]]; then
      duplex=true
      cleanup=''
    elif [[ $1 == '-c' ]]; then
      duplex=''
      cleanup=true
    elif [[ $1 == '-a' ]]; then
      duplex=true
      cleanup=true
    fi
    shift
  fi
  if [[ $# -lt 4 ]] || [[ $1 == '-h' ]]; then
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

  if ! echo "$readlen" | grep -qE '^[0-9]+$'; then
    fail "ERROR: Invalid read length \"$readlen\"."
  fi
  if ! [[ -d $outdir ]]; then
    fail "ERROR: Invalid output directory \"$outdir\"."
  fi

  print_versions $script_dir

  echo "
Parameters:
ref:     $ref
fastq1:  $fastq1
fastq2:  $fastq2
readlen: $readlen
"

  refdict=$(echo "$ref" | sed -E 's/\.fa(sta)?$//').dict
  rlenreal=$((readlen-BarcodeLen-SpacerLen))
  end_clip_start=$((rlenreal-EndClip+1))
  start=$(echo $Start | cut -d . -f 1)

  if [[ $duplex ]]; then
    if ! [[ $start ]] || [[ $start -le 56 ]]; then
      echo '===== 56 =====' && echo '===== 56 =====' >&2
      # Concatenate the 12-nt tag sequences from the paired reads and evaluate for tag quality
      $PythonCmd $script_dir/tag_to_header.py --infile1 $fastq1 --infile2 $fastq2 \
        --outfile1 $outdir/read_1.fq.smi --outfile2 $outdir/read_2.fq.smi \
        --barcode_length $BarcodeLen --spacer_length $SpacerLen
    fi
    if ! [[ $start ]] || [[ $start -le 57 ]]; then
      # Allow specifying 57.2 so only the second half is executed.
      if ! [[ $Start ]] || [[ $Start != 57.2 ]]; then
        echo '===== 57 =====' && echo '===== 57 =====' >&2
        # Align each read to the reference genome
        $BwaCmd aln $ref $outdir/read_1.fq.smi > $outdir/read_1.aln
      fi
      echo '===== 57.2 =====' && echo '===== 57.2 =====' >&2
      $BwaCmd aln $ref $outdir/read_2.fq.smi > $outdir/read_2.aln
    fi
    if ! [[ $start ]] || [[ $start -le 58 ]]; then
      echo '===== 58 =====' && echo '===== 58 =====' >&2
      # Make a single paired-end .sam file
      $BwaCmd sampe -s $ref $outdir/read_1.aln $outdir/read_2.aln \
        $outdir/read_1.fq.smi $outdir/read_2.fq.smi > $outdir/PE_reads.sam
    fi
    if ! [[ $start ]] || [[ $start -le 59 ]]; then
      echo '===== 59 =====' && echo '===== 59 =====' >&2
      # Convert to .bam format and sort by position
      $SamtoolsCmd view -Sbu $outdir/PE_reads.sam | $SamtoolsCmd sort - $outdir/PE_reads.sort
    fi
    if ! [[ $start ]] || [[ $start -le 60 ]]; then
      echo '===== 60 =====' && echo '===== 60 =====' >&2
      # Run the Python program 'ConsensusMaker.py' to collapse PCR duplicates into SSCS
      $PythonCmd $script_dir/ConsensusMaker.py --tagfile $outdir/PE_reads.tagcounts \
        --infile $outdir/PE_reads.sort.bam --outfile $outdir/SSCS.bam --min 3 --max 1000 \
        --cutoff 0.7 --Ncutoff 0.3 --readlength $rlenreal --read_type dpm --filt osn
    fi
    if ! [[ $start ]] || [[ $start -le 61 ]]; then
      echo '===== 61 =====' && echo '===== 61 =====' >&2
      # Sort the SSCS reads
      $SamtoolsCmd view -bu $outdir/SSCS.bam | $SamtoolsCmd sort - $outdir/SSCS.sort
    fi
    if ! [[ $start ]] || [[ $start -le 62 ]]; then
      echo '===== 62 =====' && echo '===== 62 =====' >&2
      # Construct DCSs from SSCSs using 'DuplexMaker.py'
      $PythonCmd $script_dir/DuplexMaker.py --infile $outdir/SSCS.sort.bam \
        --outfile $outdir/DCS_data.bam --Ncutoff 0.3 --readlength $rlenreal
    fi
  fi

  if [[ $cleanup ]]; then
    if ! [[ $start ]] || [[ $start -le 63 ]]; then
      if ! [[ -s $outdir/DCS_data.r1.fq ]] || ! [[ -s $outdir/DCS_data.r2.fq ]]; then
        fail "ERROR: Duplex sequences missing! ($outdir/DCS_data.r[12].fq)"
      fi
      if ! [[ $Start ]] || [[ $Start != 63.2 ]]; then
        echo '===== 63 =====' && echo '===== 63 =====' >&2
        # Align each DCS .fastq to the reference genome
        $BwaCmd aln $ref $outdir/DCS_data.r1.fq > $outdir/DCS_data.r1.aln
      fi
      echo '===== 63.2 =====' && echo '===== 63.2 =====' >&2
      $BwaCmd aln $ref $outdir/DCS_data.r2.fq > $outdir/DCS_data.r2.aln
    fi
    if ! [[ $start ]] || [[ $start -le 64 ]]; then
      echo '===== 64 =====' && echo '===== 64 =====' >&2
      # Make a paired-end .sam file for the DCS data
      $BwaCmd sampe -s $ref $outdir/DCS_data.r1.aln $outdir/DCS_data.r2.aln \
        $outdir/DCS_data.r1.fq $outdir/DCS_data.r2.fq > $outdir/DCS_PE.aln.sam
    fi
    if ! [[ $start ]] || [[ $start -le 65 ]]; then
      echo '===== 65 =====' && echo '===== 65 =====' >&2
      # Convert to .bam format and sort by position
      $SamtoolsCmd view -Sbu $outdir/DCS_PE.aln.sam | $SamtoolsCmd sort - $outdir/DCS_PE.aln.sort
    fi
    if ! [[ $start ]] || [[ $start -le 66 ]]; then
      echo '===== 66 =====' && echo '===== 66 =====' >&2
      # Index the final sorted DCS .bam file
      $SamtoolsCmd index $outdir/DCS_PE.aln.sort.bam
    fi
    if ! [[ $start ]] || [[ $start -le 67 ]]; then
      echo '===== 67 =====' && echo '===== 67 =====' >&2
      # Filter out unmapped reads from the final DCS .bam file
      $SamtoolsCmd view -F 4 -b $outdir/DCS_PE.aln.sort.bam > $outdir/DCS_PE.filt.bam
    fi
    if ! [[ $start ]] || [[ $start -le 68 ]]; then
      echo '===== 68 =====' && echo '===== 68 =====' >&2
      # Add readgroups field to the header of the final DCS .bam file with Picard to allow for
      # compatibility with the GATK using Picard tools.
      $JavaCmd -jar -Xmx2g $PicardDir/AddOrReplaceReadGroups.jar INPUT=$outdir/DCS_PE.filt.bam \
        OUTPUT=$outdir/DCS_PE.filt.readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default
    fi
    if ! [[ $start ]] || [[ $start -le 69 ]]; then
      if ! [[ $Start ]] || [[ $Start != 69.2 ]]; then
        echo '===== 69 =====' && echo '===== 69 =====' >&2
        # Index the final sorted DCS .bam file
        $SamtoolsCmd index $outdir/DCS_PE.filt.readgroups.bam
      fi
      if ! [[ -s $refdict ]]; then
        echo '===== 69.2 =====' && echo '===== 69.2 =====' >&2
        # rm -f $refdict
        $JavaCmd -jar $PicardDir/CreateSequenceDictionary.jar REFERENCE=$ref OUTPUT=$refdict
      fi
    fi
    if ! [[ $start ]] || [[ $start -le 70 ]]; then
      if ! [[ $Start ]] || [[ $Start != 70.2 ]]; then
        echo '===== 70 =====' && echo '===== 70 =====' >&2
        # Perform local re-alignment of the reads using GATK. First identify the genome targets for
        # local re-alignment.
        $JavaCmd -Xmx2g -jar $GatkDir/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref \
          -I $outdir/DCS_PE.filt.readgroups.bam -o $outdir/DCS_PE.filt.readgroups.intervals
      fi
      echo '===== 70.2 =====' && echo '===== 70.2 =====' >&2
      # This command is followed by the actual local re-alignment.
      $JavaCmd -Xmx2g -jar $GatkDir/GenomeAnalysisTK.jar -T IndelRealigner -R $ref \
        -I $outdir/DCS_PE.filt.readgroups.bam -o $outdir/DCS_PE.filt.readgroups.realign.bam \
        -targetIntervals $outdir/DCS_PE.filt.readgroups.intervals
    fi
    if ! [[ $start ]] || [[ $start -le 71 ]]; then
      echo '===== 71 =====' && echo '===== 71 =====' >&2
      # Perform end-trimming of DCS reads. An example command that trims five bases from both the 3'
      # and 5' ends of each read is provided.
      $JavaCmd -Xmx2g -jar $GatkDir/GenomeAnalysisTK.jar -T ClipReads \
        -I $outdir/DCS_PE.filt.readgroups.realign.bam -o $outdir/DCS-final.bam -R $ref \
        --cyclesToTrim "1-$StartClip,$end_clip_start-$rlenreal" --clipRepresentation SOFTCLIP_BASES
    fi
  fi
  echo '===== DONE =====' && echo '===== DONE =====' >&2
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
  # pipeline
  echo -en 'pipeline:\te0897da\t\t'
  if ! [[ -d $script_dir ]]; then
    echo 'MISSING'
  elif ! which git >/dev/null 2>/dev/null; then
    echo 'ERROR 1'
  else
    unset CDPATH
    cd $script_dir
    if ! git log >/dev/null 2>/dev/null; then
      echo 'ERROR 2'
    else
      git log --oneline -n 1 | grep --color=never -Eo '^\S+'
    fi
    cd - >/dev/null
  fi
  # Python
  echo -en 'Python:\t\t2.7\t\t'
  if which $PythonCmd >/dev/null 2>/dev/null; then
    $PythonCmd --version 2>&1 | sed -E 's/python\s//I'
  else
    echo 'MISSING'
  fi
  # BWA
  echo -en 'BWA:\t\t0.6.2\t\t'
  if which $BwaCmd >/dev/null 2>/dev/null; then
    $BwaCmd 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  else
    echo 'MISSING'
  fi
  # Samtools
  echo -en 'Samtools:\t0.1.18\t\t'
  if which $SamtoolsCmd >/dev/null 2>/dev/null; then
    $SamtoolsCmd 2>&1 | sed -En 's/^.*version.*\s([0-9].*)$/\1/Ip'
  else
    echo 'MISSING'
  fi
  # PySAM
  echo -en 'PySAM:\t\t0.7.5\t\t'
  if $PythonCmd -c 'import pysam' 2>/dev/null; then
    $PythonCmd -c 'import pysam; print pysam.__version__'
  elif which $PythonCmd >/dev/null 2>/dev/null; then
    echo 'MISSING'
  else
    echo 'ERROR 1'
  fi
  # BioPython
  echo -en 'BioPython:\t1.62\t\t'
  if $PythonCmd -c 'import Bio' 2>/dev/null; then
    $PythonCmd -c 'import Bio; print Bio.__version__'
  elif which $PythonCmd >/dev/null 2>/dev/null; then
    echo 'MISSING'
  else
    echo 'ERROR 1'
  fi
  if ! which $JavaCmd 2>/dev/null >/dev/null; then
    echo "ERROR: Java command \"$JavaCmd\" not found." >&2
    return
  fi
  # Picard
  echo -en 'Picard:\t\t1.107\t\t'
  if [[ -f $PicardDir/picard.jar ]]; then
    $JavaCmd -jar $PicardDir/picard.jar AddOrReplaceReadGroups --version 2>&1 >/dev/null | sed -E 's/\(.*\)//'
  elif ! [[ -f $PicardDir/AddOrReplaceReadGroups.jar ]]; then
    echo 'MISSING'
  elif [[ $($JavaCmd -jar $PicardDir/AddOrReplaceReadGroups.jar 2>&1 >/dev/null | sed -En 's/^Version:?\s//ip') ]]; then
    $JavaCmd -jar $PicardDir/AddOrReplaceReadGroups.jar 2>&1 >/dev/null | sed -En 's/^Version:?\s//ip'
  else
    echo 'ERROR 1'
  fi
  # GATK
  echo -en 'GATK:\t\t2.4-9\t\t'
  if ! [[ -f $GatkDir/GenomeAnalysisTK.jar ]]; then
    echo 'MISSING'
  else
    set +e
    version=$($JavaCmd -jar $GatkDir/GenomeAnalysisTK.jar --version 2>/dev/null)
    exit_code=$?
    set -e
    if [[ $exit_code == 0 ]]; then
      echo $version
    else
      version=$($JavaCmd -jar $GatkDir/GenomeAnalysisTK.jar 2>&1 >/dev/null | sed -En 's/^.*version\s([0-9.-]+[0-9.]).*$/\1/p')
      if [[ $version ]]; then
        echo $version
      else
        echo 'ERROR 1'
      fi
    fi
  fi
}


function fail {
  echo "$@" >&2
  exit 1
}


main "$@"
