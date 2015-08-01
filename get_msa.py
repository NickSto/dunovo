#!/usr/bin/env python
from __future__ import division
import os
import sys
import argparse
import tempfile
import subprocess
import consensus
import seqtools

OPT_DEFAULTS = {'format':'plain', 'qual':20, 'qual_format':'sanger'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('seqs', metavar='sequence', nargs='*',
    help='The alignment.')
  parser.add_argument('-i', '--input',
    help='Provide the sequences in this input file instead of as command-line arguments. '
         'Give "-" to use stdin.')
  parser.add_argument('-f', '--format', choices=('plain', 'duplex'),
    help='Input format. "plain" is a simple list of the sequences, one on each line. "duplex" is '
         'the 8-column format of the family-sorted read data from the duplex pipeline. It must be '
         'the read pairs from a single alpha/beta barcode combination (both the alpha-beta and '
         'beta-alpha strands). If "duplex" is given, you must also specify which of the four '
         'possible alignments to output with --mate and --order.')
  parser.add_argument('-m', '--mate', type=int, choices=(1, 2))
  parser.add_argument('-o', '--order', choices=('ab', 'ba'))
  parser.add_argument('-F', '--qual-format', choices=('sanger',))
  parser.add_argument('-q', '--qual', type=int,
    help='Quality threshold: Default: %(default)s')

  args = parser.parse_args(argv[1:])

  qual_thres = ' '
  if args.qual_format == 'sanger':
    qual_thres = chr(args.qual + 33)
  else:
    fail('Error: Unsupported FASTQ quality format "{}".'.format(args.qual_format))
  # Check arguments.
  if not (args.seqs or args.input):
    fail('Error: You must provide sequences either in a file with --input or as arguments.')
  elif args.seqs and args.input:
    fail('Error: You cannot provide sequences in both a file and command-line arguments.')
  if args.format == 'duplex' and not (args.mate and args.order):
    fail('Error: If the --format is duplex, you must specify a --mate and --order.')

  # Read input.
  quals = []
  if args.input:
    if args.format == 'plain':
      if args.input == '-':
        seqs = [line.strip() for line in sys.stdin]
      else:
        with open(args.input) as infile:
          seqs = [line.strip() for line in infile]
    elif args.format == 'duplex':
      if args.input == '-':
        (seqs, quals) = parse_duplex(sys.stdin, args.mate, args.order)
      else:
        with open(args.input) as infile:
          (seqs, quals) = parse_duplex(infile, args.mate, args.order)
  else:
    seqs = args.seqs

  align = make_msa(seqs)
  if quals:
    quals = seqtools.transfer_gaps_multi(quals, align, gap_char_out=' ')
  cons = consensus.get_consensus(align, quals, qual_thres=qual_thres, gapped=True)

  output = format_alignment(cons, align, quals, qual_thres=ord(qual_thres))

  for seq in output:
    print seq


def parse_duplex(infile, mate, order):
  seqs = []
  quals = []
  for line in infile:
    (bar, this_order, name1, seq1, qual1, name2, seq2, qual2) = line.rstrip('\r\n').split('\t')
    if this_order == order:
      if mate == 1:
        seqs.append(seq1)
        quals.append(qual1)
      elif mate == 2:
        seqs.append(seq2)
        quals.append(qual2)
  return seqs, quals


def make_msa(seqs):
  """Perform a multiple sequence alignment on a set of sequences.
  Uses MAFFT."""
  i = 0
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='msa.') as family_file:
    for seq in seqs:
      i+=1
      header = '>{}\n'.format(i)
      family_file.write(header)
      family_file.write(seq+'\n')
  with open(os.devnull, 'w') as devnull:
    try:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      output = subprocess.check_output(command, stderr=devnull)
    except (OSError, subprocess.CalledProcessError):
      return None
  os.remove(family_file.name)
  return read_fasta(output)


def read_fasta(fasta):
  """Quick and dirty FASTA parser. Return only the list of sequences (no names).
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  sequence = ''
  for line in fasta.splitlines():
    if line.startswith('>'):
      if sequence:
        sequences.append(sequence)
      sequence = ''
      continue
    sequence += line.strip()
  if sequence:
    sequences.append(sequence)
  return sequences


def format_alignment(cons, seqs, quals=(), qual_thres=32, id_char='.'):
  output = [cons.upper()]
  for i, seq in enumerate(seqs):
    outseq = ''
    for j, seq_base in enumerate(seq.upper()):
      if quals and seq_base != '-' and ord(quals[i][j]) < qual_thres:
        outseq += ' '
      elif cons[j] == seq_base:
        outseq += id_char
      else:
        outseq += seq_base
    output.append(outseq)
  return output


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
