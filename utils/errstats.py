#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse

ARG_DEFAULTS = {'input':sys.stdin, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """Tally statistics on errors in reads, compared to the rest of their (single-
stranded) families.
Output columns without --all-repeats:
1. barcode
2. order
3. mate
4. number of reads
5. number of errors that occurred more than once
6-end. number of errors in each read.
With --all-repeats:
1. barcode
2. order
3. mate
4. number of reads
5-end. count of occurrence of each error that occurred more than once."""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='families.msa.tsv', nargs='?', type=argparse.FileType('r'),
    help='')
  parser.add_argument('-H', '--human', action='store_true')
  parser.add_argument('-r', '--all-repeats', action='store_true')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  run(args.input, **vars(args))


def run(infile, human=False, all_repeats=False, *nargs, **kwargs):
  for family in parse(infile):
    for order in ('ab', 'ba'):
      for mate in (0, 1):
        alignment = family[order][mate]
        if alignment:
          consensus, errors, repeat_errors, new_align = compare(alignment, all_repeats=all_repeats,
                                                                make_new_align=human)
          if human:
            for seq, seq_errors in zip(new_align, errors):
              print('{} errors: {}'.format(seq, seq_errors))
            if all_repeats:
              print('{} errors: {}, repeat errors: {}\n'.format(consensus, sum(errors),
                                                                ', '.join(map(str, repeat_errors))))
            else:
              print('{} errors: {}, repeat errors: {}\n'.format(consensus, sum(errors), repeat_errors))
          elif all_repeats:
            print(family['bar'], order, mate, len(alignment), *repeat_errors, sep='\t')
          else:
            print(family['bar'], order, mate, len(alignment), repeat_errors, *errors, sep='\t')


def parse(infile):
  last_barcode = None
  family = {'bar':None, 'ab':([], []), 'ba':([], [])}
  for line in infile:
    fields = line.split('\t')
    barcode = fields[0]
    order = fields[1]
    mate = int(fields[2])-1
    seq = fields[4]
    if barcode != last_barcode:
      if last_barcode is not None:
        yield family
      family = {'bar':barcode, 'ab':([], []), 'ba':([], [])}
      last_barcode = barcode
    family[order][mate].append(seq)
  yield family


def compare(alignment, all_repeats=False, make_new_align=False):
  num_seqs = len(alignment)
  consensus = ''
  new_align = [''] * num_seqs
  errors = [0] * num_seqs
  if all_repeats:
    repeat_errors = []
  else:
    repeat_errors = 0
  for bases in zip(*alignment):
    # Tally how many of each base there are at this position, and find the "winner" (the consensus).
    votes = {'A':0, 'C':0, 'G':0, 'T':0, '-':0, 'N':0}
    cons = 'N'
    max_vote = 0
    for base in bases:
      vote = votes[base]
      vote += 1
      if vote == max_vote:
        cons = 'N'
      elif vote > max_vote:
        max_vote = vote
        cons = base
      votes[base] = vote
    # Make new alignment without consensus bases.
    for i, base in enumerate(bases):
      if base == cons:
        new_align[i] += '.'
      else:
        new_align[i] += base
        errors[i] += 1
    # How often do we see the same error more than once?
    for base, vote in votes.items():
      if base != cons and vote > 1:
        if all_repeats:
          repeat_errors.append(vote)
        else:
          repeat_errors += 1
    consensus += cons
    i += 1
  return consensus, errors, repeat_errors, new_align


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
