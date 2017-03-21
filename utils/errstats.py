#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse

ARG_DEFAULTS = {'infile':sys.stdin, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('infile', metavar='families.msa.tsv', nargs='?', type=argparse.FileType('r'),
    help='')
  parser.add_argument('-H', '--human', action='store_true')
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

  run(args.infile)


def run(infile):

  for family in parse(infile):
    for order in ('ab', 'ba'):
      for mate in (0, 1):
        alignment = family[order][mate]
        if alignment:
          consensus, errors, errors_per_seq, new_align = compare(alignment)
          for seq, errors_in_seq in zip(new_align, errors_per_seq):
            print('{} errors: {}'.format(seq, errors_in_seq))
          print('{} errors: {}'.format(consensus, errors))
          print()


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


def compare(alignment):
  errors = 0
  i = 0
  end = False
  consensus = ''
  new_align = [''] * len(alignment)
  errors_per_seq = [0] * len(alignment)
  while not end:
    votes = {'A':0, 'C':0, 'G':0, 'T':0, '-':0, 'N':0}
    for seq in alignment:
      try:
        base = seq[i]
      except IndexError:
        end = True
        break
      votes[base] += 1
    if end:
      break
    # Tally the votes.
    cons = 'N'
    max_vote = 0
    for base, vote in votes.items():
      if vote == max_vote:
        cons = 'N'
      elif vote > max_vote:
        max_vote = vote
        cons = base
    # Count errors.
    if cons != 'N':
      for base, vote in votes.items():
        if base != cons:
          errors += vote
    # Make new alignment without consensus bases.
    for j, seq in enumerate(alignment):
      if seq[i] == cons:
        new_align[j] += '.'
      else:
        new_align[j] += seq[i]
        errors_per_seq[j] += 1
    consensus += cons
    i += 1
  return consensus, errors, errors_per_seq, new_align


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
