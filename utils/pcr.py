#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse

ARG_DEFAULTS = {'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('-x', type=int, required=True,
    help='Number of reads with the error.')
  parser.add_argument('-n', type=int, required=True,
    help='Total number of reads.')
  parser.add_argument('-k', type=int, required=True,
    help='Number of PCR cycles.')
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

  print(calc(args.n, args.x, args.k))


def calc(n, x, k):
  """Calculate the equation:
  $\frac{\sum_{i=1}^k 2^i \frac{n!}{(n - x)! x!} \frac{1}{2^i}^x (1 - \frac{1}{2^i})^{n - x} }
  {\sum_{y=1}^{n-1} \sum_{i=1}^k 2^i \frac{n!}{(n-y)!y!} \frac{1}{2^i}^y (1 - \frac{1}{2^i})^{n-y}}$
  """
  numerator = summation(equation1, 1, k, n, x)
  denominator = 0
  for y in range(1, n):
    denominator += summation(equation1, 1, k, n, y)
  return numerator/denominator


def summation(function, start, end, *args):
  sum = 0
  for i in range(start, end+1):
    sum += function(i, *args)
  return sum


def equation1(i, n, x):
  mult1 = 2**i
  mult2 = factorial(n)/(factorial(n-x)*factorial(x))
  mult3 = (1/2**i)**x
  mult4 = (1-(1/2**i))**(n-x)
  return mult1 * mult2 * mult3 * mult4


def factorial(n):
  if n <= 1:
    return n
  else:
    return n * factorial(n-1)


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
