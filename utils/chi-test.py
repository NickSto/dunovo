#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import math
import errno
import logging
import argparse
import collections
import scipy.stats

ARG_DEFAULTS = {'log_file':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('expected', metavar='expected.tsv', type=argparse.FileType('r'),
    help='')
  parser.add_argument('observed', metavar='observed.tsv', type=argparse.FileType('r'),
    help='')
  parser.add_argument('-l', '--log', action='store_true')
  parser.add_argument('-L', '--log-file', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log_file, level=args.volume, format='%(message)s')
  tone_down_logger()

  dists = read_expected(args.expected)
  counts = read_observed(args.observed)
  totals = sum_counts(counts)
  all_obs_freqs = calc_freqs(counts, totals)

  for famsize, exp_freqs in dists.items():
    obs_freqs = all_obs_freqs[famsize]
    assert len(exp_freqs) == len(obs_freqs), famsize
    if len(exp_freqs) <= 1:
      continue
    # for exp, obs in zip(exp_freqs, obs_freqs):
    #   print('{:8.5f}\t{:8.5f}'.format(exp, obs))
    chi = scipy.stats.chisquare(obs_freqs, exp_freqs)
    if args.log:
      chi_stat = math.log(chi.statistic, 10)
      chi_p = math.log(chi.pvalue, 10)
    else:
      chi_stat = chi.statistic
      chi_p = chi.pvalue
    print('{}\t{}\t{:8.5f}\t{:8.5f}'.format(famsize, totals[famsize], chi_stat, chi_p))


def read_expected(expected_file):
  dists = collections.defaultdict(list)
  for line in expected_file:
    fields = line.rstrip('\r\n').split('\t')
    famsize = int(fields[1])
    prob = float(fields[3])
    dists[famsize].append(prob)
  return dists


def read_observed(observed_file):
  counts = collections.defaultdict(lambda: collections.defaultdict(int))
  for line in observed_file:
    fields = line.rstrip('\r\n').split('\t')
    famsize = int(fields[0])
    errors = int(fields[1])
    counts[famsize][errors] += 1
  return counts


def sum_counts(counts):
  totals = {}
  for famsize, observations in counts.items():
    totals[famsize] = sum(observations.values())
  return totals


def calc_freqs(counts, totals):
  freqs = collections.defaultdict(list)
  for famsize, observations in counts.items():
    total = totals[famsize]
    for num_errors in range(1, famsize//2+1):
      freqs[famsize].append(observations[num_errors] / total)
  return freqs


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
