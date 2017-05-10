#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import tempfile
import argparse
import subprocess
import collections
import multiprocessing
import distutils.spawn
from lib import simplewrap
from lib import version
from ET import phone
import seqtools

#TODO: Warn if it looks like the two input FASTQ files are the same (i.e. the _1 file was given
#      twice). Can tell by whether the alpha and beta (first and last 12bp) portions of the barcodes
#      are always identical. This would be a good thing to warn about, since it's an easy mistake
#      to make, but it's not obvious that it happened. The pipeline won't fail, but will just
#      produce pretty weird results.

REQUIRED_COMMANDS = ['mafft']
OPT_DEFAULTS = {'processes':1}
DESCRIPTION = """Read in sorted FASTQ data and do multiple sequence alignments of each family."""


def main(argv):

  wrapper = simplewrap.Wrapper()
  wrap = wrapper.wrap

  parser = argparse.ArgumentParser(description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawTextHelpFormatter)
  parser.set_defaults(**OPT_DEFAULTS)

  wrapper.width = wrapper.width - 24
  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help=wrap('The input reads, sorted into families. One line per read pair, 8 tab-delimited '
              'columns:\n'
              '1. canonical barcode\n'
              '2. barcode order ("ab" for alpha+beta, "ba" for beta-alpha)\n'
              '3. read 1 name\n'
              '4. read 1 sequence\n'
              '5. read 1 quality scores\n'
              '6. read 2 name\n'
              '7. read 2 sequence\n'
              '8. read 2 quality scores'))
  parser.add_argument('-p', '--processes', type=int,
    help=wrap('Number of worker subprocesses to use. Must be at least 1. Default: %(default)s.'))
  parser.add_argument('--phone-home', action='store_true',
    help='Report helpful usage data to the developer, to better understand the use cases and '
         'performance of the tool. The only data which will be recorded is the name and version of '
         'the tool, the size of the input data, the time taken to process it, and the IP address '
         'of the machine running it. No parameters or filenames are sent. All the reporting and '
         'recording code is available at https://github.com/NickSto/ET.')
  parser.add_argument('--test', action='store_true',
    help='If reporting usage data, mark this as a test run.')
  parser.add_argument('-v', '--version', action='version', version=str(version.get_version()),
    help=wrap('Print the version number and exit.'))

  args = parser.parse_args(argv[1:])

  start_time = time.time()
  if args.phone_home:
    run_id = phone.send_start(__file__, version.get_version(), test=args.test)

  assert args.processes > 0, '-p must be greater than zero'

  # Check for required commands.
  missing_commands = []
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      missing_commands.append(command)
  if missing_commands:
    fail('Error: Missing commands: "'+'", "'.join(missing_commands)+'".')

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  # Open all the worker processes.
  workers = open_workers(args.processes)

  # Main loop.
  # This processes whole duplexes (pairs of strands) at a time for a future option to align the
  # whole duplex at a time.
  stats = {'duplexes':0, 'time':0, 'pairs':0, 'runs':0, 'aligned_pairs':0}
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
        output, run_stats, current_worker_i = delegate(workers, stats, duplex, barcode)
        process_results(output, run_stats, stats)
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
  output, run_stats, current_worker_i = delegate(workers, stats, duplex, barcode)
  process_results(output, run_stats, stats)

  # Do one last loop through the workers, reading the remaining results and stopping them.
  # Start at the worker after the last one processed by the previous loop.
  start = current_worker_i + 1
  for i in range(len(workers)):
    worker_i = (start + i) % args.processes
    worker = workers[worker_i]
    output, run_stats = worker.recv()
    process_results(output, run_stats, stats)
    worker.send(None)

  if infile is not sys.stdin:
    infile.close()

  end_time = time.time()
  run_time = int(end_time - start_time)

  # Final stats on the run.
  sys.stderr.write('Processed {pairs} read pairs in {duplexes} duplexes.\n'.format(**stats))
  if stats['aligned_pairs'] > 0:
    per_pair = stats['time'] / stats['aligned_pairs']
    per_run = stats['time'] / stats['runs']
    sys.stderr.write('{:0.3f}s per pair, {:0.3f}s per run.\n'.format(per_pair, per_run))
    sys.stderr.write('{}s total time.\n'.format(run_time))

  if args.phone_home:
    stats['align_time'] = stats['time']
    del stats['time']
    phone.send_end(__file__, version.get_version(), run_id, run_time, stats, test=args.test)


def open_workers(num_workers):
  """Open the required number of worker processes."""
  workers = []
  for i in range(num_workers):
    parent_pipe, child_pipe = multiprocessing.Pipe()
    process = multiprocessing.Process(target=worker_function, args=(child_pipe,))
    process.start()
    workers.append(parent_pipe)
  return workers


def worker_function(child_pipe):
  while True:
    args = child_pipe.recv()
    if args is None:
      break
    child_pipe.send(process_duplex(*args))


def delegate(workers, stats, duplex, barcode):
  worker_i = stats['duplexes'] % len(workers)
  worker = workers[worker_i]
  if stats['duplexes'] >= len(workers):
    output, run_stats = worker.recv()
  else:
    output, run_stats = '', {}
  stats['duplexes'] += 1
  args = (duplex, barcode)
  worker.send(args)
  return output, run_stats, worker_i


def process_duplex(duplex, barcode):
  output = ''
  run_stats = {'time':0, 'runs':0, 'aligned_pairs':0}
  orders = duplex.keys()
  if len(duplex) == 0 or None in duplex:
    return '', {}
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
    start = time.time()
    alignment = align_family(family, mate)
    # Compile statistics.
    elapsed = time.time() - start
    pairs = len(family)
    #logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
    if pairs > 1:
      run_stats['time'] += elapsed
      run_stats['runs'] += 1
      run_stats['aligned_pairs'] += pairs
    if alignment is None:
      pass  #logging.warning('Error aligning family {}/{} (read {}).'.format(barcode, order, mate))
    else:
      output += format_msa(alignment, barcode, order, mate)
  return output, run_stats


def align_family(family, mate):
  """Do a multiple sequence alignment of the reads in a family and their quality scores."""
  mate = str(mate)
  assert mate == '1' or mate == '2'
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


def process_results(output, run_stats, stats):
  """Process the outcome of a duplex run.
  Print the aligned output and sum the stats from the run with the running totals."""
  for key, value in run_stats.items():
    stats[key] += value
  if output:
    sys.stdout.write(output)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
