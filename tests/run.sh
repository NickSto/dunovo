#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)

USAGE="Usage: \$ $(basename $0) [options] [test1 [test2]]"


function main {

  do_all=true
  verbose=true
  # Run the requested tests
  for arg in "$@"; do
    # Check for options
    #TODO: option to keep test data at end instead of removing it.
    if [[ ${arg:0:1} == '-' ]]; then
      case "$arg" in
        -h)
          echo "$USAGE" >&2
          echo "Currently valid tests:" >&2
          list_tests >&2
          exit 1;;
        -q)
          verbose='';;
        -v)
          verbose=true;;
        *)
          echo "Unrecognized option \"$arg\"." >&2;;
      esac
      continue
    fi
    # Execute valid tests (if they're existing functions).
    if [[ $(type -t $arg) == function ]]; then
      do_all=''
      if [[ $verbose ]]; then
        $arg
      else
        $arg 2>/dev/null
      fi
    else
      echo "Unrecognized test \"$arg\"." >&2
      do_all=''
    fi
  done

  # If no tests were specified in arguments, do all tests.
  if [[ $do_all ]]; then
    if [[ $verbose ]]; then
      all
    else
      all 2>/dev/null
    fi
  fi
}

function fail {
  echo "$@" >&2
  exit 1
}

function list_tests {
  while read declare f test; do
    # Filter out functions that aren't tests.
    if echo "$initial_declarations" | grep -qF 'declare -f '"$test"; then
      continue
    else
      echo "$test"
    fi
  done < <(declare -F)
}

# Capture a list of all functions defined before the tests, to tell which are actual functions
# and which are tests.
initial_declarations=$(declare -F)

########## Functional tests ##########

# Do all tests.
function all {
  sscs
  sscs_r1
  sscs_p3
}

# sscs.py defaults
function sscs {
  echo -e "\tsscs.py ::: families.in.tsv:"
  python "$dirname/../sscs.py" "$dirname/families.in.tsv" | diff -s - "$dirname/families.out.tsv"
}

# sscs.py no minimum number of reads per family
function sscs_r1 {
  echo -e "\tsscs.py -r 1 ::: families.in.tsv:"
  python "$dirname/../sscs.py" -r 1 "$dirname/families.in.tsv" | diff -s - "$dirname/families-r1.out.tsv"
}

# sscs.py parallelized
function sscs_p3 {
  echo -e "\tsscs.py -p 3 ::: families.in.tsv:"
  python "$dirname/../sscs.py" -p 3 "$dirname/families.in.tsv" | paste - - | sort \
    | diff -s - <(cat "$dirname/families.out.tsv" | paste - - | sort)
}

main "$@"
