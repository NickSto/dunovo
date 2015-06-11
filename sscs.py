#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse

OPT_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build single-strand consensus sequences from read families.
Pipe sorted reads into stdin."""

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')

  args = parser.parse_args(argv[1:])

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  current_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 7:
      continue
    (barcode, name1, seq1, qual1, name2, seq2, qual2) = fields
    if barcode != current_barcode:
      current_barcode = barcode
      if family:
        process_family(family, barcode)
      family = []
    family.append((name1, seq1, qual1, name2, seq2, qual2))

  if infile is not sys.stdin:
    infile.close()


def process_family(family, barcode):
  print '>'+barcode
  for (name1, seq1, qual1, name2, seq2, qual2) in family:
    print '\t'+name1
    print '\t'+name2


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
