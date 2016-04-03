#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import logging
import argparse
import tempfile
import subprocess
import multiprocessing
import consensus
import swalign

ARG_DEFAULTS = {'bar_len':24, 'win_len':5, 'shift':3, 'processes':1, 'loglevel':logging.ERROR}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Try to match barcodes with sequencing errors.
Match based on a small window in the middle of each half of the barcode.
Then it will align all the unique barcodes which match and then print the similarity of each to the
consensus."""
EPILOG = """This will print each kmer observed, the barcodes which contained it, and their
similarities. The output is 4 tab-delimited columns: 1. whether the kmer was in the first or second
half of the barcode (0 for first half, 1 for second) 2. the kmer 3. the barcode 4. its similarity to
the consensus"""

# Algorithm from Paul Medvedev (email from 2015-12-16)

def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('infile', metavar='families.tsv', nargs='?',
    help='Input file (sorted output of make-barcodes.awk).')
  parser.add_argument('-n', '--num-barcodes', type=int,
    help='Only read in this many different barcodes.')
  parser.add_argument('-c', '--consensus', action='store_true',
    help='Include consensus sequences in the output. They will appear the same as normal barcodes, '
         'but they will be printed before each set of barcodes matching a kmer. (So you can filter '
         'them out by looking for when either column 1 or 2 change, then discard the line after '
         'the change.')
  parser.add_argument('-b', '--bar-len', type=int,
    help='Barcode length. Default: %(default)s')
  parser.add_argument('-w', '--win-len', type=int,
    help='Window (k-mer) size. Default: %(default)s')
  parser.add_argument('-s', '--shift', type=int,
    help='Bases to shift the window (number of k-mers to check). Default: %(default)s')
  parser.add_argument('-q', '--quiet', dest='loglevel', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='loglevel', action='store_const', const=logging.INFO)
  parser.add_argument('--debug', dest='loglevel', action='store_const', const=logging.DEBUG)
  parser.add_argument('-p', '--processes', type=int,
    help='Number of worker processes to use. Default: %(default)s')

  args = parser.parse_args(argv[1:])

  assert args.processes > 0, '-p must be greater than zero'
  logging.basicConfig(stream=sys.stderr, level=args.loglevel, format='%(message)s')

  starts = calc_starts(args.bar_len, args.win_len, args.shift)

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  logging.info('Beginning to read in data.')
  # For each window sequence (kmer), build a set of barcodes which contained it, in any of the shift
  # positions. Do this for both halves of the barcode (independently, at the moment).
  kmer_dicts = [{}, {}]
  last_barcode = None
  barcode_count = 0
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 8:
      logging.warn('Line contains incorrect number of fields.')
      continue
    barcode = fields[0]
    # Only do it for each unique barcode (in the sorted output, there will be runs of lines with
    # the same barcode).
    if barcode == last_barcode:
      continue
    barcode_count += 1
    # for each half of the barcode
    for kmer_dict, start in zip(kmer_dicts, starts):
      # for each shift position (trying kmers at each of args.shift offsets)
      for i in range(args.shift):
        kmer = barcode[start+i:start+i+args.win_len]
        kmer_set = kmer_dict.get(kmer, set())
        kmer_set.add(barcode)
        kmer_dict[kmer] = kmer_set
    last_barcode = barcode
    if args.num_barcodes and barcode_count >= args.num_barcodes:
      break

  if infile is not sys.stdin:
    infile.close()

  workers = open_workers(args.processes)

  # Analyze the groups of barcodes that contained each kmer:
  # Multiple sequence align all the barcodes in a each, call a consensus, then smith-waterman
  # align each barcode to that consensus to measure their similarity to it.
  run_num = 0
  for dict_num, kmer_dict in enumerate(kmer_dicts):
    # Each half of the barcode (one dict per).
    for kmer, barcodes_set in kmer_dict.items():
      # Each set of barcodes which share a kmer.
      barcodes = list(barcodes_set)
      results = delegate(workers, run_num, dict_num, kmer, barcodes)
      if results:
        process_results(*results, print_consensus=args.consensus)
      run_num += 1

  # Do one last loop through the workers, reading the remaining results and stopping them.
  # Start at the worker after the last one processed by the previous loop.
  start_i = (run_num + 1) % len(workers)
  for i in range(len(workers)):
    worker_i = (start_i + i) % args.processes
    worker = workers[worker_i]
    results = worker.recv()
    if results:
      process_results(*results, print_consensus=args.consensus)
    worker.send(None)


def calc_starts(bar_len, win_len, shift):
  half_len = bar_len//2
  assert win_len < half_len, 'Window length must be less than half the barcode length.'
  # Place the window right in the middle of the first half of the barcode.
  # Offset is where it should start.
  offset = (half_len-win_len)/2
  # Move it backward by half the shift length so that the average kmer start is at the offset
  # calculated above.
  start1 = int(offset - shift/2)
  start2 = start1 + half_len
  return start1, start2


def process_results(dict_num, kmer, consensus_seq, barcodes, similarities, print_consensus=False):
  if print_consensus:
    print(dict_num, kmer, consensus_seq, 1.0, sep='\t')
  for barcode, similarity in zip(barcodes, similarities):
    print(dict_num, kmer, barcode, similarity, sep='\t')


def open_workers(num_workers):
  """Open the required number of worker processes."""
  workers = []
  for i in range(num_workers):
    parent_pipe, child_pipe = multiprocessing.Pipe()
    process = multiprocessing.Process(target=worker_function, args=(child_pipe,))
    process.start()
    workers.append(parent_pipe)
  return workers


def delegate(workers, run_num, dict_num, kmer, barcodes):
  worker_i = run_num % len(workers)
  worker = workers[worker_i]
  if run_num >= len(workers):
    logging.info('Parent: Trying to receive results from worker..')
    results = worker.recv()
  else:
    results = None
  args = (dict_num, kmer, barcodes)
  logging.info('Parent: Sending new data to worker..')
  worker.send(args)
  return results


##### HAPPENS IN CHILD PROCESSES #####

def worker_function(child_pipe):
  while True:
    # logging.info('Worker: Listening for new data from parent..')
    args = child_pipe.recv()
    if args is None:
      break
    # logging.info('Worker: Sending results back to parent..')
    child_pipe.send(process_barcodes(*args))


def process_barcodes(dict_num, kmer, barcodes):
  """Perform a multiple sequence alignment on a set of barcodes and parse the result.
  Uses MAFFT."""
  # If there's only one barcode, we don't have to do an alignment.
  if len(barcodes) == 1:
    return dict_num, kmer, barcodes[0], barcodes, [1.0]
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.') as family_file:
    for i, barcode in enumerate(barcodes):
      family_file.write('>{}\n'.format(i))
      family_file.write(barcode+'\n')
  with open(os.devnull, 'w') as devnull:
    try:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      output = subprocess.check_output(command, stderr=devnull)
    except (OSError, subprocess.CalledProcessError):
      return None
  os.remove(family_file.name)
  alignment = read_fasta(output, upper=True)
  consensus_seq = consensus.get_consensus(alignment)
  similarities = []
  for barcode in barcodes:
    similarities.append(get_similarity(consensus_seq, barcode))
  return dict_num, kmer, consensus_seq, barcodes, similarities


def read_fasta(fasta, upper=False):
  """Quick and dirty FASTA parser. Return a list of the sequences only (no names)."""
  sequences = []
  sequence = ''
  for line in fasta.splitlines():
    if line.startswith('>'):
      if upper:
        sequence = sequence.upper()
      if sequence:
        sequences.append(sequence)
      sequence = ''
      continue
    sequence += line.strip()
  if upper:
    sequence = sequence.upper()
  if sequence:
    sequences.append(sequence)
  return sequences


def get_similarity(seq1, seq2):
  align = swalign.smith_waterman(seq1, seq2)
  logging.debug(align.target+'\n'+align.query)
  return align.matches / len(align.query)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
