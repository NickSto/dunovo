#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import subprocess

TMP_DIRNAME = 'tmp-sscs'
OPT_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build single-strand consensus sequences from read families. Pipe sorted reads into
stdin. Prints single-strand consensus sequences to stdout. The sequence names are BARCODE.MATE, e.g.
"CTCAGATAACATACCTTATATGCA.1"."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')

  args = parser.parse_args(argv[1:])

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  family_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 7:
      continue
    (barcode, name1, seq1, qual1, name2, seq2, qual2) = fields
    if barcode != family_barcode:
      if family:
        process_family(family, family_barcode)
      family_barcode = barcode
      family = []
    family.append((name1, seq1, qual1, name2, seq2, qual2))
  if family:
    process_family(family, barcode)

  if infile is not sys.stdin:
    infile.close()


def process_family(family, barcode):
  if not os.path.isdir(TMP_DIRNAME):
    os.mkdir(TMP_DIRNAME)
  if len(family) == 1:
    (name1, seq1, qual1, name2, seq2, qual2) = family[0]
    print '>'+barcode+'.1'
    print seq1
    print '>'+barcode+'.2'
    print seq2
    return
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
      subprocess.call(['mafft', family_path], stdout=align_file, stderr=devnull)
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
