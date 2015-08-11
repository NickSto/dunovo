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
import distutils.spawn
import seqtools

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

  # Open all the worker processes, if we're using more than one.
  workers = None
  if args.processes > 1:
    workers = open_workers(args.processes, args)

  # Main loop.
  # This processes whole duplexes (pairs of strands) at a time for a future option to align the
  # whole duplex at a time.
  stats = {'time':0, 'pairs':0, 'runs':0, 'families':0}
  all_pairs = 0
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
        process_duplex(duplex, barcode, workers=workers, processes=args.processes, stats=stats)
        duplex = collections.OrderedDict()
      barcode = this_barcode
      order = this_order
      family = []
    pair = {'name1': name1, 'seq1':seq1, 'qual1':qual1, 'name2':name2, 'seq2':seq2, 'qual2':qual2}
    family.append(pair)
    all_pairs += 1
  # Process the last family.
  duplex[order] = family
  # sys.stderr.write('processing {}: {} orders ({}) [last]\n'.format(barcode, len(duplex),
  #                  '/'.join([str(len(duplex[order])) for order in duplex])))
  process_duplex(duplex, barcode, workers=workers, processes=args.processes, stats=stats)

  if args.processes > 1:
    close_workers(workers)
    compile_results(workers)
    delete_tempfiles(workers)

  if infile is not sys.stdin:
    infile.close()

  if not args.stats_file:
    return

  # Final stats on the run.
  logging.info('Processed {} read pairs and {} multi-pair families.'
               .format(all_pairs, stats['runs']))
  per_pair = stats['time'] / stats['pairs']
  per_run = stats['time'] / stats['runs']
  logging.info('{:0.3f}s per pair, {:0.3f}s per run.'.format(per_pair, per_run))


def open_workers(num_workers, args):
  """Open the required number of worker processes."""
  script_path = os.path.realpath(sys.argv[0])
  workers = []
  for i in range(num_workers):
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
    outfile = tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.out.part.')
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=outfile)
    worker = {'proc':process, 'outfile':outfile, 'stats':stats_subfile}
    workers.append(worker)
  return workers


def delegate(worker, duplex, barcode):
  """Send a family to a worker process."""
  for order, family in duplex.items():
    for pair in family:
      line = '{}\t{}\t{name1}\t{seq1}\t{qual1}\t{name2}\t{seq2}\t{qual2}\n'.format(barcode, order,
                                                                                   **pair)
      worker['proc'].stdin.write(line)


def close_workers(workers):
  for worker in workers:
    worker['outfile'].close()
    worker['proc'].stdin.close()


def compile_results(workers):
  for worker in workers:
    worker['proc'].wait()
    with open(worker['outfile'].name, 'r') as outfile:
      for line in outfile:
        sys.stdout.write(line)


def delete_tempfiles(workers):
  for worker in workers:
    os.remove(worker['outfile'].name)
    if worker['stats']:
      os.remove(worker['stats'])


def process_duplex(duplex, barcode, workers=None, processes=1, stats=None):
  # Are we the controller process or a worker?
  stats['families'] += 1
  if processes > 1:
    i = stats['families'] % len(workers)
    worker = workers[i]
    delegate(worker, duplex, barcode)
    return
  # We're a worker. Actually process the family.
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
      logging.warning('Error aligning family {}/{} (read {}).'.format(barcode, order, mate))
    else:
      print_msa(alignment, barcode, order, mate)


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
  logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
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
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as family_file:
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


def print_msa(align, barcode, order, mate, outfile=sys.stdout):
  for sequence in align:
    outfile.write('{bar}\t{order}\t{mate}\t{name}\t{seq}\t{qual}\n'
                  .format(bar=barcode, order=order, mate=mate, **sequence))


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
