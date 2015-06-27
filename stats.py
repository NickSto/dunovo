#!/usr/bin/env python
from __future__ import division
import os
import sys
import ctypes
import argparse
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
align = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'align.so'))

OPT_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.msa.tsv', nargs='?',
    help='The --msa output of sscs.py. Will read from stdin if not provided.')
  parser.add_argument('-s', '--str',
    help='default: %(default)s')

  args = parser.parse_args(argv[1:])

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  consensus = None
  barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 3:
      continue
    (this_barcode, name, seq) = fields
    if fields[1] == 'CONSENSUS':
      if family and consensus:
        process_family(barcode, consensus, family)
      barcode = this_barcode
      consensus = seq
      family = []
    else:
      family.append(seq)
  if family and consensus:
    process_family(barcode, consensus, family)

  if infile is not sys.stdin:
    infile.close()


#TODO: Maybe print the number of N's in the consensus?
def process_family(barcode, consensus, family):
  diffs = get_diffs(consensus, family)
  for (i, diff) in enumerate(diffs):
    print "{}\t{}\t{}".format(barcode, diff, family[i].upper())


def get_diffs(consensus, family):
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    c_family[i] = ctypes.c_char_p(seq)
  align.get_diffs_simple.restype = ctypes.POINTER(ctypes.c_int * len(c_family))
  diffs = align.get_diffs_simple(c_consensus, c_family, len(c_family))
  return diffs.contents


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
