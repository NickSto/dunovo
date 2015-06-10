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

  parser.add_argument('positional1', metavar='dispname',
    help='')
  parser.add_argument('-s', '--str',
    help='default: %(default)s')

  args = parser.parse_args(argv[1:])


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
