#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import re
import sys
import errno
import logging
import argparse

ARG_DEFAULTS = {'families':sys.stdin, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """Filter out reads from the families.tsv file based on barcode. This will print the
input to stdout, omitting lines with barcodes that fail the filter(s)."""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('families', metavar='families.tsv', type=argparse.FileType('r'),
    help='Any tab-delimited file where the first column is the barcode (both alpha and beta). '
         'For instance, the output of make-families.awk or align_families.py.')
  parser.add_argument('-r', '--repeats', type=int,
    help='Longest allowed single-base repeat. Default: No limit.')
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

  run(args.families, sys.stdout, args.repeats)


def run(families, output, max_repeats):
  for line in families:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) <= 1:
      logging.warn('Invalid line ({} fields)'.format(len(fields)))
      continue
    barcode = fields[0]
    halflen = len(barcode)//2
    halves = [barcode[halflen:], barcode[:halflen]]
    # Apply filters.
    passes = True
    longest_longest_repeat = 0
    if max_repeats is not None:
      bases = set(barcode)
      regexes = make_regexes(bases)
      for half in halves:
        longest_repeat = find_longest_repeat(half, regexes)
        longest_longest_repeat = max(longest_longest_repeat, longest_repeat)
        if longest_repeat > max_repeats:
          passes = False
    if passes:
      output.write(line)
    else:
      logging.info('Filtered out {} (repeat {} > {})'
                   .format(barcode, longest_longest_repeat, max_repeats))


def make_regexes(bases):
  regexes = []
  for base in bases:
    regexes.append(re.compile(base+'+'))
  return regexes


def find_longest_repeat(barcode, regexes):
  longest_repeat = 0
  for regex in regexes:
    repeats = regex.findall(barcode)
    for repeat in repeats:
      longest_repeat = max(longest_repeat, len(repeat))
  return longest_repeat


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
