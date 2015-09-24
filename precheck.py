#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import getreads

OPT_DEFAULTS = {'tag_len':12, 'const_len':5}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Print statistics on the raw duplex sequencing reads."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infiles', nargs='*',
    help='The input files. Give the first and second mate files if FASTQ, e.g. "reads_1.fq '
         'reads_2.fq". If SAM, give one file. Omit to read from stdin.')
  parser.add_argument('-f', '--format', choices=('fastq',))
  parser.add_argument('-t', '--tag-length', dest='tag_len', type=int)
  parser.add_argument('-c', '--constant-length', dest='const_len', type=int)

  args = parser.parse_args(argv[1:])

  filetype = args.format
  if len(args.infiles) == 0:
    if not filetype:
      fail('Error: must provide a --format if no input files are given (reading from stdin).')
  if len(args.infiles) == 1:
    if not filetype:
      if args.infiles[0].endswith('.sam'):
        filetype = 'sam'
      else:
        fail('Error: if giving 1 input file, it must end in .sam or you need to specify a '
             '--format.')
    if filetype != 'sam':
      fail('Error: valid formats for 1 input file: sam.')
  elif len(args.infiles) == 2:
    if not filetype:
      if ((args.infiles[0].endswith('.fq') or args.infiles[0].endswith('.fastq')) and
          (args.infiles[1].endswith('.fq') or args.infiles[1].endswith('.fastq'))):
        filetype = 'fastq'
      else:
        fail('Error: if giving 2 input files, they must both end in .fq or .fastq or you need to '
             'specify a --format.')
    if filetype != 'fastq':
      fail('Error: valid formats for 2 input files: fastq.')

  infileh1 = None
  infileh2 = None
  if len(args.infiles) == 0:
    infileh1 = sys.stdin
  else:
    infileh1 = open(args.infiles[0])
  reader1 = getreads.getparser(infileh1, filetype=filetype).parser()
  if len(args.infiles) == 2:
    infileh2 = open(args.infiles[1])
    reader2 = getreads.getparser(infileh2, filetype=filetype).parser()

  total_reads = 0
  barcodes = {}
  while True:
    if filetype == 'fastq':
      try:
        read1 = reader1.next()
        read2 = reader2.next()
      except StopIteration:
        break
    barcode = read1.seq[:args.tag_len] + read2.seq[:args.tag_len]
    if barcode not in barcodes:
      barcodes[barcode] = 1
    else:
      barcodes[barcode] += 1
    total_reads += 1

  if infileh1 and infileh1 is not sys.stdin:
    infileh1.close()
  if infileh2 and infileh2 is not sys.stdin:
    infileh2.close()

  print 'Total read pairs:\t'+str(total_reads)
  print 'Unique barcodes:\t'+str(len(barcodes))
  print 'Read pairs per barcode:\t'+str(total_reads/len(barcodes))



def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
