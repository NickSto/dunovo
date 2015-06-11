#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import argparse
import subprocess
import distutils.spawn

TMP_DIRNAME = 'tmp-sscs'
REQUIRED_COMMANDS = ('mafft', 'em_cons')
OPT_DEFAULTS = {'processes':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build single-strand consensus sequences from read families. Pipe sorted reads into
stdin. Prints single-strand consensus sequences to stdout. The sequence names are BARCODE.MATE, e.g.
"CTCAGATAACATACCTTATATGCA.1"."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')
  parser.add_argument('-s', '--stats-file',
    help='Print statistics on the run to this file. Use "-" to print to stderr.')
  # parser.add_argument('-p', '--processes', type=int,
  #   help='Number of processes to use. If > 1, launches this many worker subprocesses. '
  #        'Default: %(default)s.')

  args = parser.parse_args(argv[1:])

  # Check for required commands.
  missing_commands = []
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      missing_commands.append(command)
  if missing_commands:
    fail('Error: Missing commands "'+'", "'.join(missing_commands)+'".')

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  if args.stats_file:
    if args.stats_file == '-':
      logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(message)s')
    else:
      logging.basicConfig(filename=args.stats_file, filemode='w', level=logging.INFO,
                          format='%(message)s')
  else:
    logging.disable(logging.CRITICAL)

  total_time = 0
  total_pairs = 0
  total_runs = 0
  all_pairs = 0
  family = []
  family_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 7:
      continue
    (barcode, name1, seq1, qual1, name2, seq2, qual2) = fields
    # If the barcode has changed, we're in a new family.
    # Process the reads we've previously gathered as one family and start a new family.
    if barcode != family_barcode:
      if family:
        (elapsed, pairs) = process_family(family, family_barcode)
        if pairs > 1:
          total_time += elapsed
          total_pairs += pairs
          total_runs += 1
      family_barcode = barcode
      family = []
    family.append((name1, seq1, qual1, name2, seq2, qual2))
    all_pairs += 1
  # Process the last family.
  if family:
    (elapsed, pairs) = process_family(family, family_barcode)
    if pairs > 1:
      total_time += elapsed
      total_pairs += pairs
      total_runs += 1

  # Final stats on the run.
  logging.info('Processed {} read pairs and {} multi-pair families.'.format(all_pairs, total_runs))
  per_pair = total_time / total_pairs
  per_run = total_time / total_runs
  logging.info('{:0.2f}s per pair, {:0.2f}s per run.'.format(per_pair, per_run))

  if infile is not sys.stdin:
    infile.close()


def process_family(family, barcode):
  if not os.path.isdir(TMP_DIRNAME):
    os.mkdir(TMP_DIRNAME)
  start = time.time()
  pairs = len(family)
  if pairs == 1:
    (name1, seq1, qual1, name2, seq2, qual2) = family[0]
    print '>'+barcode+'.1'
    print seq1
    print '>'+barcode+'.2'
    print seq2
  else:
    align_path = make_msa(family, 1)#, base=barcode+'.1')
    consensus = get_consensus(align_path)#, base=barcode+'.1')
    if consensus is not None:
      print '>'+barcode+'.1'
      print consensus
    align_path = make_msa(family, 2)#, base=barcode+'.2')
    consensus = get_consensus(align_path)#, base=barcode+'.2')
    if consensus is not None:
      print '>'+barcode+'.2'
      print consensus
  end = time.time()
  elapsed = end - start
  logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
  return (elapsed, pairs)


def make_msa(family, mate, base='family'):
  """Perform a multiple sequence alignment on a set of sequences.
  Uses MAFFT."""
  family_path = os.path.join(TMP_DIRNAME, base+'.fa')
  with open(family_path, 'w') as family_file:
    for pair in family:
      if mate == 1:
        name = pair[0]
        seq = pair[1]
      else:
        name = pair[3]
        seq = pair[4]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  align_path = os.path.join(TMP_DIRNAME, base+'.align.fa')
  with open(align_path, 'w') as align_file:
    with open(os.devnull, 'w') as devnull:
      subprocess.call(['mafft', '--nuc', '--quiet', family_path], stdout=align_file, stderr=devnull)
  return align_path


def get_consensus(align_path, base='family'):
  """Make a consensus from a multiple sequence alignment file and return the
  consensus sequence as a string.
  Uses the EMBOSS em_cons command."""
  # Note on em_cons output:
  # It may always be lowercase, but maybe not. It can contain "N", and possibly "-".
  cons_path = os.path.join(TMP_DIRNAME, base+'.cons.fa')
  with open(os.devnull, 'w') as devnull:
    subprocess.call(['em_cons', '-sequence', align_path, '-outseq', cons_path], stderr=devnull)
  if not os.path.exists(cons_path):
    return None
  return read_fasta(cons_path)


def read_fasta(fasta_path):
  """Read a FASTA file, return the sequence.
  Uses a very narrow definition of FASTA: That returned by the "em_cons" command."""
  seq_lines = []
  at_header = True
  with open(fasta_path) as fasta_file:
    for line in fasta_file:
      if at_header and line.startswith('>'):
        at_header = False
        continue
      seq_lines.append(line.strip())
  return "".join(seq_lines)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
