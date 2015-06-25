#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import distutils.spawn

REQUIRED_COMMANDS = ['mafft', 'em_cons']
OPT_DEFAULTS = {'min_reads':3, 'processes':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build single-strand consensus sequences from read families. Pipe sorted reads into
stdin. Prints single-strand consensus sequences in FASTA to stdout. The sequence names are
BARCODE.MATE, e.g. "CTCAGATAACATACCTTATATGCA.1"."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families. One line per read pair, 7 tab-delimited columns: '
         '1. barcode (concatenated alpha+beta), 2. read 1 name, 3. read 1 sequence, 4. read 1 '
         'quality scores, 5. read 2 name, 6. read 2 sequence, 7. read 2 quality scores.')
  parser.add_argument('-r', '--min-reads', type=int,
    help='The minimum number of reads required to form a family. Families with fewer reads will '
         'be skipped. Default: %(default)s.')
  parser.add_argument('-m', '--msa', action='store_true',
    help='Print the multiple sequence alignment as well as the consensus. Instead of printing '
         'FASTA, it will print a tab-delimited format with 3 columns: 1. barcode, 2. read name, 3. '
         'aligned read sequence. One line will be the consensus sequence, named "CONSENSUS". One '
         'alignment and consensus will be printed for each read in the pair.')
  parser.add_argument('-s', '--stats-file',
    help='Print statistics on the run to this file. Use "-" to print to stderr.')
  parser.add_argument('-p', '--processes', type=int,
    help='Number of processes to use. If > 1, launches this many worker subprocesses. '
         'Default: %(default)s.')
  parser.add_argument('-S', '--slurm', action='store_true',
    help='If -p > 1, prepend sub-commands with "srun -C new".')

  args = parser.parse_args(argv[1:])

  assert args.processes > 0

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

  stats = {'time':0, 'pairs':0, 'runs':0, 'families':0}
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
      process_family(family, family_barcode, args, workers=workers, stats=stats)
      family_barcode = barcode
      family = []
    family.append((name1, seq1, qual1, name2, seq2, qual2))
    all_pairs += 1
  # Process the last family.
  process_family(family, family_barcode, args, workers=workers, stats=stats)

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
    arguments = gather_args(sys.argv, args.infile)
    command.extend(arguments)
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


def gather_args(args, infile, excluded_flags=('-S', '--slurm'),
                excluded_args=('-p', '--processes', '-s', '--stats-file')):
  """Take the full list of command-line arguments and return only the ones which
  should be passed to worker processes.
  Excludes the 0th argument (the command name), the input filename ("infile"), all
  arguments in "excluded_flags", and all arguments in "excluded_args" plus the
  arguments which follow."""
  out_args = []
  skip = True
  for arg in args:
    if skip:
      skip = False
      continue
    if arg in excluded_flags:
      continue
    if arg in excluded_args:
      skip = True
      continue
    if arg == infile:
      continue
    out_args.append(arg)
  return out_args


def delegate(worker, family, barcode):
  """Send a family to a worker process."""
  for pair in family:
    line = barcode+'\t'+'\t'.join(pair)+'\n'
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


def process_family(family, barcode, args, workers=None, stats=None):
  # Pass if family doesn't contain minimum # of reads.
  if len(family) < args.min_reads:
    return
  stats['families'] += 1
  # Are we the controller process or a worker?
  if args.processes > 1:
    i = stats['families'] % len(workers)
    worker = workers[i]
    delegate(worker, family, barcode)
    return
  # We're a worker. Actually process the family.
  start = time.time()
  pairs = len(family)
  if pairs == 1:
    (name1, seq1, qual1, name2, seq2, qual2) = family[0]
    if args.msa:
      print '{bar}\tCONSENSUS\t{seq}\n{bar}\t{name}\t{seq}'.format(bar=barcode, name=name1, seq=seq1)
      print '{bar}\tCONSENSUS\t{seq}\n{bar}\t{name}\t{seq}'.format(bar=barcode, name=name2, seq=seq2)
    else:
      print '>'+barcode+'.1'
      print seq1
      print '>'+barcode+'.2'
      print seq2
  else:
    output = 'consensus'
    if args.msa:
      output = 'msa'
    align_and_cons(family, barcode, 1, output=output)
    align_and_cons(family, barcode, 2, output=output)
  end = time.time()
  elapsed = end - start
  logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
  if stats and pairs > 1:
    stats['time'] += elapsed
    stats['pairs'] += pairs
    stats['runs'] += 1
  return


def align_and_cons(family, barcode, mate, output='consensus'):
  """Do the bioinformatic work: make the multiple sequence alignment and consensus
  sequence.
  "mate" is "1" or "2" and determines which read in the pair to process.
  "output" is "consensus" to print the consensus sequence in FASTA format, or
    "msa" to print the full multiple sequence alignment in tab-delimited format."""
  align_path = make_msa(family, mate)
  if output == 'msa':
    msa = read_fasta(align_path)
  consensus = get_consensus(align_path)
  if output == 'consensus' and consensus is not None:
    print '>{}.{}'.format(barcode, mate)
    print consensus
  elif output == 'msa':
    print '{}\t{}\t{}'.format(barcode, 'CONSENSUS', consensus)
    for sequence in msa:
      print '{}\t{}\t{}'.format(barcode, sequence['name'], sequence['seq'])


def make_msa(family, mate):
  """Perform a multiple sequence alignment on a set of sequences.
  Uses MAFFT."""
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as family_file:
    for pair in family:
      if mate == 1:
        name = pair[0]
        seq = pair[1]
      else:
        name = pair[3]
        seq = pair[4]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as align_file:
    with open(os.devnull, 'w') as devnull:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      subprocess.call(command, stdout=align_file, stderr=devnull)
  os.remove(family_file.name)
  return align_file.name


def get_consensus(align_path):
  """Make a consensus from a multiple sequence alignment file and return the
  consensus sequence as a string.
  Uses the EMBOSS em_cons command."""
  # Note on em_cons output:
  # It may always be lowercase, but maybe not. It can contain "N", and possibly "-".
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as cons_file:
    cons_path = cons_file.name
  with open(os.devnull, 'w') as devnull:
    command = ['em_cons', '-sequence', align_path, '-outseq', cons_path]
    subprocess.call(command, stderr=devnull)
  os.remove(align_path)
  if os.path.getsize(cons_path) == 0:
    os.remove(cons_path)
    return None
  else:
    seqs = read_fasta(cons_path)
    os.remove(cons_path)
    if seqs:
      return seqs[0]['seq']


def read_fasta(fasta_path):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences. Each is a dict of 'name' and 'seq'."""
  sequences = []
  seq_lines = []
  seq_name = None
  with open(fasta_path) as fasta_file:
    for line in fasta_file:
      if line.startswith('>'):
        if seq_lines:
          sequences.append({'name':seq_name, 'seq':''.join(seq_lines)})
        seq_lines = []
        seq_name = line.rstrip('\r\n')[1:]
        continue
      seq_lines.append(line.strip())
  if seq_lines:
    sequences.append({'name':seq_name, 'seq':''.join(seq_lines)})
  return sequences


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
