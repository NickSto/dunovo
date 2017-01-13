#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse
import subprocess

ARG_DEFAULTS = {'initial_barcodes':20, 'mapq_thres':25, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('initial_barcodes', metavar='barcodes to try', type=int, nargs='?',
    help='')
  parser.add_argument('-m', '--mapq', type=int)
  parser.add_argument('-r', '--random', action='store_true')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  logging.info('Reading random barcodes from border-families.txt..')
  # cat border-families.txt | paste - - | shuf --random-source=border-families.txt | head -n $initial_barcodes
  if args.random:
    shuf_cmd = ('shuf',)
  else:
    shuf_cmd = ('shuf', '--random-source=border-families.txt')
  process = make_pipeline(('cat', 'border-families.txt'),
                          ('paste', '-', '-'),
                          shuf_cmd,
                          ('head', '-n', str(args.initial_barcodes)))
  families_by_barcode = {}
  for line in process.stdout:
    fields = line.split()
    family = {}
    family['count1'] = int(fields[0])
    family['count2'] = int(fields[3])
    barcode = fields[1]
    family['barcode'] = barcode
    families_by_barcode[barcode] = family

  logging.info('Reading barcodes.fq to find read names..')
  hits = 0
  families_by_read_name = {}
  line_num = 0
  with open('barcodes.fq', 'rU') as barcodes_fq:
    for line in barcodes_fq:
      line_num += 1
      if line_num % 4 == 1:
        read_name = line[1:].rstrip('\r\n')
      elif line_num % 4 == 2:
        seq = line.rstrip('\r\n')
        family = families_by_barcode.get(seq)
        if family:
          hits += 1
          family['read_name'] = read_name
          families_by_read_name[read_name] = family
  logging.info('hits: {}'.format(hits))

  logging.info('Reading barcodes.bam to find similar barcodes..')
  hits = 0
  neighbors_by_read_name = {}
  # samtools view -f 256 barcodes.bam | awk '$1 == '$read_name' && $5 > 25 {print $3}'
  process = subprocess.Popen(('samtools', 'view', '-f', '256', 'barcodes.bam'), stdout=subprocess.PIPE)
  for line in process.stdout:
    fields = line.split()
    mapq = int(fields[4])
    if mapq >= args.mapq_thres:
      read_name = fields[0]
      family = families_by_read_name.get(read_name)
      if family:
        hits += 1
        read_name2 = fields[2]
        neighbor = {'read_name':read_name2}
        neighbors = family.get('neighbors', [])
        neighbors.append(neighbor)
        family['neighbors'] = neighbors
        neighbors_by_read_name[read_name2] = neighbor
  logging.info('hits: {}'.format(hits))

  logging.info('Reading barcodes.fq to find sequences of similar barcodes..')
  hits = 0
  line_num = 0
  with open('barcodes.fq', 'rU') as barcodes_fq:
    for line in barcodes_fq:
      line_num += 1
      if line_num % 4 == 1:
        read_name = line[1:].rstrip('\r\n')
        neighbor = neighbors_by_read_name.get(read_name)
      if line_num % 4 == 2 and neighbor:
        seq = line.rstrip('\r\n')
        neighbor['barcode'] = seq
  logging.info('hits: {}'.format(hits))

  logging.info('Printing results..')
  for family in families_by_barcode.values():
    neighbors = family.get('neighbors')
    if not neighbors:
      continue
    print('{barcode}\t{count1}\t{count2}\t{read_name}'.format(**family))
    for neighbor in neighbors:
      print(neighbor['barcode'])


def make_pipeline(*commands):
  processes = []
  for command in commands:
    if not processes:
      processes.append(subprocess.Popen(command, stdout=subprocess.PIPE))
    else:
      processes.append(subprocess.Popen(command, stdin=processes[-1].stdout, stdout=subprocess.PIPE))
  processes[0].stdout.close()
  return processes[-1]


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
