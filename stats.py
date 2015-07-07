#!/usr/bin/env python
from __future__ import division
import os
import sys
import ctypes
import argparse
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
align = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'align.so'))

BINS = 10

OPT_DEFAULTS = {'bins':10}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


#TODO: Report diffs as a percentage of sequence length (or bin length), since they can vary.
def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('stat', choices=('diffs', 'diffs-binned', 'seqlen'),
    help='The type of statistics to compute and print.')
  parser.add_argument('infile', metavar='read-families.msa.tsv', nargs='?',
    help='The --msa output of sscs.py. Will read from stdin if not provided.')
  parser.add_argument('-b', '--bins', type=int,
    help='The number of bins to segment reads into when doing "diffs-binned".')

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
        process_family(args, barcode, consensus, family)
      barcode = this_barcode
      consensus = seq
      family = []
    else:
      family.append(seq)
  if family and consensus:
    process_family(args, barcode, consensus, family)

  if infile is not sys.stdin:
    infile.close()


#TODO: Maybe print the number of N's in the consensus?
def process_family(args, barcode, consensus, family):
  if args.stat == 'diffs':
    diffs = get_diffs(consensus, family)
    for (i, diff) in enumerate(diffs):
      print "{}\t{}\t{}".format(barcode, diff, family[i].upper())
  elif args.stat == 'diffs-binned':
    diffs = get_diffs_binned(consensus, family, args.bins)
    if diffs is None:
      return
    for (i, bins) in enumerate(diffs):
      sys.stdout.write(barcode+'\t')
      for diff in bins.contents:
        sys.stdout.write(str(diff)+'\t')
      sys.stdout.write(family[i].upper()+'\n')
  elif args.stat == 'seqlen':
    seq_lens = get_seq_lens(family)
    print barcode+'\t'+'\t'.join(map(str, seq_lens))


def get_seq_lens(family):
  seq_lens = []
  for seq in family:
    seq_lens.append(len(seq))
  return seq_lens


def get_diffs(consensus, family):
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    c_family[i] = ctypes.c_char_p(seq)
  align.get_diffs_simple.restype = ctypes.POINTER(ctypes.c_int * len(c_family))
  diffs = align.get_diffs_simple(c_consensus, c_family, len(c_family))
  return diffs.contents


def get_diffs_binned(consensus, family, bins):
  seq_len = None
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    if seq_len:
      if seq_len != len(seq):
        return None
    else:
      seq_len = len(seq)
    c_family[i] = ctypes.c_char_p(seq)
  int_array_pointer = ctypes.POINTER(ctypes.c_int * bins)
  align.get_diffs_binned.restype = ctypes.POINTER(int_array_pointer * len(c_family))
  diffs = align.get_diffs_binned(c_consensus, c_family, len(c_family), seq_len, bins)
  return diffs.contents


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
