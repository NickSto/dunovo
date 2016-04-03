#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import numpy
import bisect
import random
import shutil
import tempfile
import argparse
import subprocess
import fastqreader

ARG_DEFAULTS = {'read_len':100, 'frag_len':400, 'n_frags':1000, 'seq_error':0.001, 'pcr_error':0.001,
                'cycles':25, 'indel_rate':0.15, 'extension_rate':0.3, 'seed':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""

RAW_DISTRIBUTION = (
  #  0     1     2     3     4     5     6     7     8     9
  # Low singletons, but then constant drop-off. From pML113 (see 2015-09-28 report).
  #  0,  100,   36,   31,   27,   22,   17,   12,    7,  4.3,
  #2.4,  1.2,  0.6,  0.3,  0.2, 0.15,  0.1, 0.07, 0.05, 0.03,
  # High singletons, but then a second peak around 10. From Christine plasmid (2015-10-06 report).
     0,  100, 5.24, 3.67, 3.50, 3.67, 3.85, 4.02, 4.11, 4.20,
  4.17, 4.10, 4.00, 3.85, 3.69, 3.55, 3.38, 3.15, 2.92, 2.62,
  2.27, 2.01, 1.74, 1.56, 1.38, 1.20, 1.02, 0.85,
)


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('ref', metavar='ref.fa', nargs='?',
    help='Reference sequence.')
  parser.add_argument('-r', '--read-len', type=int,
    help='Default: %(default)s')
  parser.add_argument('-f', '--frag-len', type=int,
    help='Default: %(default)s')
  parser.add_argument('-F', '--n-frags', type=int,
    help='Default: %(default)s')
  parser.add_argument('-s', '--seq-error', type=float,
    help='Sequencing error rate per base (0-1 proportion, not percent). Default: %(default)s')
  parser.add_argument('-p', '--pcr-error', type=float,
    help='PCR error rate per base (0-1 proportion, not percent). Default: %(default)s')
  parser.add_argument('-c', '--cycles', type=int,
    help='Number of PCR cycles to simulate. Default: %(default)s')
  parser.add_argument('-i', '--indel-rate', type=float,
    help='Fraction of errors which are indels. Default: %(default)s')
  parser.add_argument('-E', '--extension-rate', type=float,
    help='Probability an indel is extended. Default: %(default)s')
  parser.add_argument('--frag-file',
    help='A file of fragments already generated with wgsim. Use this instead of generating a new '
         'one.')
  parser.add_argument('-S', '--seed', type=int,
    help='Seed. Default: %(default)s')

  args = parser.parse_args(argv[1:])

  assert args.ref or args.frag_file

  # Overview: wgsim "reads" (fragments) from the reference with 0 mutation rate, then wgsim actual
  # duplex reads from each fragment.

  #TODO: Check for wgsim on the PATH.

  # Create a temporary director to do our work in. Then work inside a try so we can finally remove
  # the directory no matter what exceptions are encountered.
  tmpdir = tempfile.mkdtemp(prefix='wgdsim.')
  try:
    # Step 1: Use wgsim to create fragments from the reference.
    if args.frag_file:
      frag_file = args.frag_file
    else:
      frag_file = os.path.join(tmpdir, 'fragments.fq')
      #TODO: Check exit status
      run_command('wgsim', '-e', '0', '-r', '0', '-d', '0', '-R', args.indel_rate, '-N', args.n_frags,
                  '-X', args.extension_rate, '-1', args.frag_len, args.ref, frag_file, os.devnull)

    # NOTE: Coordinates here are 0-based (0 is the first base in the sequence).
    count = 0
    extended_dist = extend_dist(RAW_DISTRIBUTION)
    proportional_dist = compile_dist(extended_dist)
    for fragment in fastqreader.FastqReadGenerator(frag_file):
      # print('>'+fragment.name, fragment.seq, '+', fragment.qual, sep='\n')
      read1_orig = fragment.seq[:args.read_len]
      # Step 2: Determine how many reads to produce from each fragment.
      # - Use random.random() and divide the range 0-1 into segments of sizes proportional to
      #   the likelihood of each family size.
      # bisect.bisect() finds where an element belongs in a sorted list, returning the index.
      # proportional_dist is just such a sorted list, with values from 0 to 1.
      nreads = bisect.bisect(proportional_dist, random.random())

      # Step 3: Introduce PCR errors.
      # - Determine the mutations and their frequencies.
      #   - Could get frequency from the cycle of PCR it occurs in.
      #     - Important to have PCR errors shared between reads.
      # - For each read, determine which mutations it contains.
      #   - Use random.random() < mut_freq.

      # Step 4: Introduce sequencing errors.
      read1s = []
      for i in range(nreads):
        read1 = read1_orig
        for mutation in generate_mutations(args.read_len, args.seq_error, args.indel_rate,
                                           args.extension_rate):
          #TODO: Note that an indel earlier in the read will mean later indels won't be inserted at
          #      the originally intended coordinate. Fix?
          read1 = apply_mutation(mutation, read1)
        read1s.append(read1)
      #TODO: Add barcodes and invariant sequences.

      # Print family.
      for read1 in read1s:
        count += 1
        print('>'+str(count), read1, '+', 'I' * len(read1), sep='\n')

  finally:
    shutil.rmtree(tmpdir)


def run_command(*command):
  devnull = open(os.devnull, 'w')
  try:
    exit_status = subprocess.call(map(str, command), stderr=devnull)
  except OSError:
    exit_status = None
  finally:
    devnull.close()
  return exit_status


def extend_dist(raw_dist, exponent=1.25, min_prob=0.00001, max_len_mult=2):
  """Add an exponentially decreasing tail to the distribution.
  It takes the final value in the distribution and keeps dividing it by
  "exponent", adding each new value to the end. It will not add probabilities
  smaller than "min_prob" or extend the length of the list by more than
  "max_len_mult" times."""
  extended_dist = list(raw_dist)
  final_sum = sum(raw_dist)
  value = raw_dist[-1]
  value /= exponent
  while value/final_sum >= min_prob and len(extended_dist) < len(raw_dist)*max_len_mult:
    extended_dist.append(value)
    final_sum += value
    value /= exponent
  return extended_dist


def compile_dist(raw_dist):
  """Turn the human-readable list of probabilities defined at the top into
  proportional probabilities.
  E.g. [10, 5, 5] -> [0.5, 0.75, 1.0]"""
  proportional_dist = []
  final_sum = sum(raw_dist)
  current_sum = 0
  for magnitude in raw_dist:
    current_sum += magnitude
    proportional_dist.append(current_sum/final_sum)
  return proportional_dist


#TODO: Rewrite.
"""No more simulating every mutation and whether it affected one of the final reads.
Instead, just take each of the final reads and ask whether a mutation occurred in any one of its
ancestors. But unfortunately we can't just step through each previous cycle and generate errors.
What we care about is errors shared between reads. To address that, at each cycle we have to ask
whether the read shares that ancestor with any of the other reads. Then we have to propagate each of
those errors to all shared reads. Basically, we need to work backwards, building a tree of how each
of the reads is related. Essentially, coalescent theory.
Here's the plan:
The tree will be composed of nested lists of lists, as suggested here:
https://stackoverflow.com/questions/3009935/looking-for-a-good-python-tree-data-structure
We'll go through the reads, and for each one, build a chain upward. At each cycle, we'll ask whether
this ancestor is shared with any of the other reads (Can use numpy.random.binomial(num_reads, p)
where p = 1/2**cycle (the probability of sharing an ancestor at that level). If it returns 1 or
greater, then we randomly choose a different read and join their chains.)
After building the tree, we go back down it, starting from the root, simulating mutations. Can use
recursion to propagate a mutation down both branches when they occur. Do this for every base.
"""
def make_pcr_errors(cycles, error_rate, read_len, n_reads, min_prob):
  # For each base..
  for i in range(read_len):
    # For each cycle..
    for cycle in range(1, cycles+1):
      replications = 2**cycle
      if n_reads * 1/replications < min_prob:
        continue
      n_muts = numpy.random.binomial(replications, error_rate)


def generate_mutations(read_len, error_rate, indel_rate, extension_rate):
  """Generate all the mutations that occur over the length of a read."""
  for i in range(read_len):
    if random.random() < error_rate:
      mtype, alt = make_mutation(indel_rate, extension_rate)
      yield {'coord':i, 'type':mtype, 'alt':alt}
  # Allow for an insertion after the last base.
  if random.random() < error_rate:
    mtype, alt = make_mutation(indel_rate, extension_rate)
    if mtype == 'ins':
      yield {'coord':i+1, 'type':mtype, 'alt':alt}


def make_mutation(indel_rate, extension_rate):
  """Simulate a random mutation."""
  # Is it an indel?
  if random.random() < indel_rate:
    # Is it an insertion or deletion? Decide, then initialize it.
    if random.random() < 0.5:
      mtype = 'del'
      alt = 1
    else:
      mtype = 'ins'
      alt = get_rand_base()
    # Extend the indel as long as the extension rate allows.
    while random.random() < extension_rate:
      if mtype == 'ins':
        alt += get_rand_base()
      else:
        alt += 1
  else:
    # What is the new base for the SNP?
    mtype = 'snp'
    alt = get_rand_base()
  return mtype, alt


def apply_mutation(mut, seq):
  i = mut['coord']
  if mut['type'] == 'snp':
    # Replace the base at "coord".
    new_seq = seq[:i] + mut['alt'] + seq[i+1:]
  else:
    # Indels are handled by inserting or deleting bases starting *before* the base at "coord".
    # This goes agains the VCF convention, but it allows deleting the first and last base, as well
    # as inserting before and after the sequence without as much special-casing.
    if mut['type'] == 'ins':
      new_seq = seq[:i] + mut['alt'] + seq[i:]
    else:
      new_seq = seq[:i] + seq[i+mut['alt']:]
  return new_seq


def get_rand_base():
  return random.choice(('A', 'C', 'G', 'T'))


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
