#!/usr/bin/env python
from __future__ import division
import os
import sys
import random
import argparse

OPT_DEFAULTS = {'fraction':0.1, 'seed':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')
  parser.add_argument('-f', '--fraction', type=float,
    help='Fraction of families to output. Default: %(default)s')
  parser.add_argument('-s', '--seed', type=int,
    help='Random number generator seed. Default: %(default)s')

  args = parser.parse_args(argv[1:])

  random.seed(args.seed)

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  last_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if not fields:
      continue
    barcode = fields[0]
    if barcode != last_barcode:
      if random.random() <= args.fraction:
        sys.stdout.write(''.join(family))
      family = []
    family.append(line)
    last_barcode = barcode

  if infile is not sys.stdin:
    infile.close()


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
