#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import collections
import multiprocessing
import distutils.spawn
import seqtools

#TODO: Warn if it looks like the two input FASTQ files are the same (i.e. the _1 file was given
#      twice). Can tell by whether the alpha and beta (first and last 12bp) portions of the barcodes
#      are always identical. This would be a good thing to warn about, since it's an easy mistake
#      to make, but it's not obvious that it happened. The pipeline won't fail, but will just
#      produce pretty weird results.

REQUIRED_COMMANDS = ['mafft']
OPT_DEFAULTS = {'processes':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Read in sorted FASTQ data and do multiple sequence alignments of each family."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families. One line per read pair, 8 tab-delimited columns: '
         '1. canonical barcode, 2. barcode order ("ab" for alpha+beta, "ba" for beta-alpha) 3. '
         'read 1 name, 4. read 1 sequence, 5. read 1 quality scores, 6. read 2 name, 7. read 2 '
         'sequence, 8. read 2 quality scores.')
  parser.add_argument('-s', '--stats-file',
    help='Print statistics on the run to this file. Use "-" to print to stderr.')
  parser.add_argument('-p', '--processes', type=int,
    help='Number of processes to use. If > 1, launches this many worker subprocesses. '
         'Default: %(default)s.')
  parser.add_argument('-S', '--slurm', action='store_true',
    help='If -p > 1, prepend sub-commands with "srun -C new".')

  args = parser.parse_args(argv[1:])

  assert args.processes > 0, '-p must be greater than zero'

  # Check for required commands.
  missing_commands = []
  if args.slurm:
    REQUIRED_COMMANDS.append('srun')
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      missing_commands.append(command)
  if missing_commands:
    fail('Error: Missing commands: "'+'", "'.join(missing_commands)+'".')

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

  # Open all the worker processes.
  workers = open_workers(args.processes)

  # Main loop.
  # This processes whole duplexes (pairs of strands) at a time for a future option to align the
  # whole duplex at a time.
  stats = {'duplexes':0, 'time':0, 'pairs':0, 'runs':0}
  current_worker_i = 0
  duplex = collections.OrderedDict()
  family = []
  barcode = None
  order = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 8:
      continue
    (this_barcode, this_order, name1, seq1, qual1, name2, seq2, qual2) = fields
    # If the barcode or order has changed, we're in a new family.
    # Process the reads we've previously gathered as one family and start a new family.
    if this_barcode != barcode or this_order != order:
      duplex[order] = family
      # If the barcode is different, we're at the end of the whole duplex. Process the it and start
      # a new one. If the barcode is the same, we're in the same duplex, but we've switched strands.
      if this_barcode != barcode:
        # sys.stderr.write('processing {}: {} orders ({})\n'.format(barcode, len(duplex),
        #                  '/'.join([str(len(duplex[order])) for order in duplex])))
        result, current_worker_i = delegate(workers, stats, duplex, barcode)
        sys.stdout.write(result)
        duplex = collections.OrderedDict()
      barcode = this_barcode
      order = this_order
      family = []
    pair = {'name1': name1, 'seq1':seq1, 'qual1':qual1, 'name2':name2, 'seq2':seq2, 'qual2':qual2}
    family.append(pair)
    stats['pairs'] += 1
  # Process the last family.
  duplex[order] = family
  # sys.stderr.write('processing {}: {} orders ({}) [last]\n'.format(barcode, len(duplex),
  #                  '/'.join([str(len(duplex[order])) for order in duplex])))
  result, current_worker_i = delegate(workers, stats, duplex, barcode)
  sys.stdout.write(result)

  # Do one last loop through the workers, reading the remaining results and stopping them.
  # Start at the worker after the last one processed by the previous loop.
  start = current_worker_i + 1
  for i in range(len(workers)):
    worker_i = (start + i) % args.processes
    worker = workers[worker_i]
    sys.stdout.write(worker['pipe'].recv())
    worker['pipe'].send(None)

  if infile is not sys.stdin:
    infile.close()

  if not args.stats_file:
    return

  # Final stats on the run.
  logging.info('Processed {pairs} read pairs in {duplexes} duplexes.'.format(**stats))
  per_pair = stats['time'] / stats['pairs']
  per_run = stats['time'] / stats['runs']
  logging.info('{:0.3f}s per pair, {:0.3f}s per run.'.format(per_pair, per_run))


def open_workers(num_workers):
  """Open the required number of worker processes."""
  workers = []
  for i in range(num_workers):
<<<<<<< HEAD
    if args.slurm:
      command = ['srun', '-C', 'new', 'python', script_path]
    else:
      command = ['python', script_path]
    stats_subfile = None
    if args.stats_file:
      if args.stats_file == '-':
        stats_subfile = '-'
      else:
        stats_subfile = "{}.{}.log".format(args.stats_file, i)
      command.extend(['-s', stats_subfile])
    outfile = tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.')
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=outfile)
    worker = {'proc':process, 'outfile':outfile, 'stats':stats_subfile}
=======
    worker = {}
    parent_pipe, child_pipe = multiprocessing.Pipe()
    worker['pipe'] = parent_pipe
    worker['process'] = multiprocessing.Process(target=worker_function, args=(child_pipe,))
    worker['process'].start()
>>>>>>> align_families.py: First draft of multiprocessing parallel implementation.
    workers.append(worker)
  return workers


def worker_function(pipe):
  while True:
    data = pipe.recv()
    if data is None:
      break
    args, kwargs = data
    pipe.send(process_duplex(*args, **kwargs))


def delegate(workers, stats, duplex, barcode):
  stats['duplexes'] += 1
  worker_i = stats['duplexes'] % len(workers)
  worker = workers[worker_i]
  args = (duplex, barcode)
  kwargs = {'stats':stats}
  if stats['duplexes'] >= len(workers):
    result = worker['pipe'].recv()
  else:
    result = ''
  worker['pipe'].send((args, kwargs))
  return result, worker_i


def process_duplex(duplex, barcode, stats=None):
  orders = duplex.keys()
  if len(duplex) == 0 or None in duplex:
    return
  elif len(duplex) == 1:
    # If there's only one strand in the duplex, just process the first mate, then the second.
    combos = ((1, orders[0]), (2, orders[0]))
  elif len(duplex) == 2:
    # If there's two strands, process in a criss-cross order:
    # strand1/mate1, strand2/mate2, strand1/mate2, strand2/mate1
    combos = ((1, orders[0]), (2, orders[1]), (2, orders[0]), (1, orders[1]))
  else:
    raise AssertionError('Error: More than 2 orders in duplex {}: {}'.format(barcode, orders))
  for mate, order in combos:
    family = duplex[order]
    alignment = align_family(family, mate, stats)
    if alignment is None:
      #logging.warning('Error aligning family {}/{} (read {}).'.format(barcode, order, mate))
      return ''
    else:
      return format_msa(alignment, barcode, order, mate)


def align_family(family, mate, stats):
  """Do a multiple sequence alignment of the reads in a family and their quality scores."""
  mate = str(mate)
  assert mate == '1' or mate == '2'
  start = time.time()
  # Do the multiple sequence alignment.
  seq_alignment = make_msa(family, mate)
  if seq_alignment is None:
    return None
  # Transfer the alignment to the quality scores.
  seqs = [read['seq'] for read in seq_alignment]
  quals_raw = [pair['qual'+mate] for pair in family]
  qual_alignment = seqtools.transfer_gaps_multi(quals_raw, seqs, gap_char_out=' ')
  # Package them up in the output data structure.
  alignment = []
  for aligned_seq, aligned_qual in zip(seq_alignment, qual_alignment):
    alignment.append({'name':aligned_seq['name'], 'seq':aligned_seq['seq'], 'qual':aligned_qual})
  elapsed = time.time() - start
  pairs = len(family)
  #logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
  if pairs > 1:
    stats['time'] += elapsed
    stats['pairs'] += pairs
    stats['runs'] += 1
  return alignment


def make_msa(family, mate):
  """Perform a multiple sequence alignment on a set of sequences and parse the result.
  Uses MAFFT."""
  mate = str(mate)
  assert mate == '1' or mate == '2'
  if len(family) == 0:
    return None
  elif len(family) == 1:
    # If there's only one read pair, there's no alignment to be done (and MAFFT won't accept it).
    return [{'name':family[0]['name'+mate], 'seq':family[0]['seq'+mate]}]
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.') as family_file:
    for pair in family:
      name = pair['name'+mate]
      seq = pair['seq'+mate]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  with open(os.devnull, 'w') as devnull:
    try:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      output = subprocess.check_output(command, stderr=devnull)
    except (OSError, subprocess.CalledProcessError):
      return None
  os.remove(family_file.name)
  return read_fasta(output, is_file=False, upper=True)


def read_fasta(fasta, is_file=True, upper=False):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences. Each is a dict of 'name' and 'seq'.
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  sequence = ''
  seq_name = None
  if is_file:
    with open(fasta) as fasta_file:
      fasta_lines = fasta_file.readlines()
  else:
    fasta_lines = fasta.splitlines()
  for line in fasta_lines:
    if line.startswith('>'):
      if upper:
        sequence = sequence.upper()
      if sequence:
        sequences.append({'name':seq_name, 'seq':sequence})
      sequence = ''
      seq_name = line.rstrip('\r\n')[1:]
      continue
    sequence += line.strip()
  if upper:
    sequence = sequence.upper()
  if sequence:
    sequences.append({'name':seq_name, 'seq':sequence})
  return sequences


def format_msa(align, barcode, order, mate, outfile=sys.stdout):
  output = ''
  for sequence in align:
    output += '{bar}\t{order}\t{mate}\t{name}\t{seq}\t{qual}\n'.format(bar=barcode, order=order,
                                                                       mate=mate, **sequence)
  return output


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
