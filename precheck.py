#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import getreads

OPT_DEFAULTS = {'tag_len':12, 'const_len':5, 'min_reads':3}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Print statistics on the raw duplex sequencing reads."""
EPILOG = """This seems to take about 200 bytes of memory per (12bp) tag."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile1', metavar='reads_1.fq',
    help='The first mates in the read pairs.')
  parser.add_argument('infile2', metavar='reads_2.fq',
    help='The second mates in the read pairs.')
  parser.add_argument('-t', '--tag-length', dest='tag_len', type=int)
  parser.add_argument('-c', '--constant-length', dest='const_len', type=int)
  parser.add_argument('-m', '--min-reads', type=int,
    help='The minimum number of reads required in each single-stranded family. Default: '
         '%(default)s')
  parser.add_argument('-v', '--validate', action='store_true',
    help='Check the id\'s of the reads to make sure the correct reads are mated into pairs.')

  args = parser.parse_args(argv[1:])

  with open(args.infile1) as infileh1:
    with open(args.infile2) as infileh2:
      barcodes = read_files(infileh1, infileh2, tag_len=args.tag_len, validate=args.validate)

  stats = get_stats(barcodes, tag_len=args.tag_len, min_reads=args.min_reads)
  print_stats(stats, min_reads=args.min_reads)


def read_files(infileh1, infileh2, tag_len=12, validate=False):
  reader1 = getreads.getparser(infileh1, filetype='fastq').parser()
  reader2 = getreads.getparser(infileh2, filetype='fastq').parser()
  barcodes = {}
  while True:
    try:
      read1 = reader1.next()
      read2 = reader2.next()
    except StopIteration:
      break
    if validate and read1.id != read2.id:
      raise getreads.FormatError('Read pair mismatch: "{}" and "{}"'.format(read1.id, read2.id))
    alpha = read1.seq[:tag_len]
    beta  = read2.seq[:tag_len]
    barcode = alpha + beta
    if barcode not in barcodes:
      barcodes[barcode] = 1
    else:
      barcodes[barcode] += 1
  return barcodes


def get_stats(barcodes, tag_len=12, min_reads=3):
  passed_sscs = 0
  duplexes = 0
  passed_duplexes = 0
  singletons = 0
  total_reads = 0
  for barcode, count in barcodes.items():
    total_reads += count
    if count == 1:
      singletons += 1
    if count >= min_reads:
      passed_sscs += 1
    alpha = barcode[:tag_len]
    beta  = barcode[tag_len:]
    reverse = beta + alpha
    if reverse in barcodes:
      duplexes += 1
      if count >= min_reads and barcodes[reverse] >= min_reads:
        passed_duplexes += 1
  # Each full duplex ends up being counted twice. Halve it to get the real total.
  stats = {
    'reads':total_reads,
    'barcodes':len(barcodes),
    'singletons':singletons,
    'avg_reads':total_reads/len(barcodes),
    'duplexes':duplexes//2,
    'passed_sscs':passed_sscs,
    'passed_duplexes':passed_duplexes//2,
  }
  return stats


def print_stats(stats, min_reads=3):
  print """Total read pairs:\t{reads}
Unique barcodes:\t{barcodes}
Avg # of read pairs per barcode:\t{avg_reads}
Singletons:\t{singletons}
Barcodes with reverse (other strand) present:\t{duplexes}
Passing threshold of {min_reads} reads per single-strand consensus:
\tSingle-strand consensus sequences:\t{passed_sscs}
\tDuplex consensus sequences:\t{passed_duplexes}""".format(min_reads=min_reads, **stats)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
