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
import consensus
import swalign

SANGER_START = 33
SOLEXA_START = 64
OPT_DEFAULTS = {'min_reads':3, 'processes':1, 'qual':20, 'qual_format':'sanger'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build consensus sequences from read aligned families. Prints duplex consensus
sequences in FASTA to stdout. The sequence ids are BARCODE.MATE, e.g. "CTCAGATAACATACCTTATATGCA.1",
where "BARCODE" is the input barcode, and "MATE" is "1" or "2" as an arbitrary designation of the
two reads in the pair. The id is followed by the count of the number of reads in the two families
(one from each strand) that make up the duplex, in the format READS1/READS2. If the duplex is
actually a single-strand consensus because the matching strand is missing, only one number is
listed."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The output of align_families.py. 6 columns: 1. (canonical) barcode. 2. order ("ab" or '
         '"ba"). 3. mate ("1" or "2"). 4. read name. 5. aligned sequence. 6. aligned quality '
         'scores.')
  parser.add_argument('-r', '--min-reads', type=int,
    help='The minimum number of reads (from each strand) required to form a single-strand '
         'consensus. Strands with fewer reads will be skipped. Default: %(default)s.')
  parser.add_argument('-q', '--qual', type=int,
    help='Base quality threshold. Bases below this quality will not be counted. '
         'Default: %(default)s.')
  parser.add_argument('-F', '--qual-format', choices=('sanger', 'solexa'),
    help='FASTQ quality score format. Sanger scores are assumed to begin at \'{}\' ({}). Default: '
         '%(default)s.'.format(SANGER_START, chr(SANGER_START)))
  parser.add_argument('--incl-sscs', action='store_true',
    help='When outputting duplex consensus sequences, include reads without a full duplex (missing '
         'one strand). The result will just be the single-strand consensus of the remaining read.')
  parser.add_argument('-s', '--sscs-file',
    help='Save single-strand consensus sequences in this file (FASTA format). Currently does not '
         'work when in parallel mode.')
  parser.add_argument('-l', '--log', metavar='LOG_FILE', dest='stats_file',
    help='Print statistics on the run to this file. Use "-" to print to stderr.')
  parser.add_argument('-p', '--processes', type=int,
    help='Number of processes to use. If > 1, launches this many worker subprocesses. Note: if '
         'this option is used, no output will be generated until the end of the entire run, so no '
         'streaming is possible. Default: %(default)s.')
  parser.add_argument('-S', '--slurm', action='store_true',
    help='If --processes > 1, prepend sub-commands with "srun -C new".')

  args = parser.parse_args(argv[1:])

  assert args.processes > 0, '-p must be greater than zero'
  # Make dict of process_family() parameters that don't change between families.
  static = {}
  static['processes'] = args.processes
  static['incl_sscs'] = args.incl_sscs
  static['min_reads'] = args.min_reads
  if args.sscs_file:
    static['sscs_fh'] = open(args.sscs_file, 'w')
  if args.qual_format == 'sanger':
    static['qual_thres'] = chr(args.qual + SANGER_START)
  elif args.qual_format == 'solexa':
    static['qual_thres'] = chr(args.qual + SOLEXA_START)
  else:
    fail('Error: unrecognized --qual-format.')

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

  stats = {'time':0, 'reads':0, 'runs':0, 'families':0}
  all_reads = 0
  duplex = collections.OrderedDict()
  family = []
  barcode = None
  order = None
  mate = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 6:
      continue
    (this_barcode, this_order, this_mate, name, seq, qual) = fields
    this_mate = int(this_mate)
    # If the barcode or order has changed, we're in a new single-stranded family.
    # Process the reads we've previously gathered as one family and start a new family.
    if this_barcode != barcode or this_order != order or this_mate != mate:
      duplex[(order, mate)] = family
      # We're at the end of the duplex pair if the barcode changes or if the order changes without
      # the mate changing, or vice versa (the second read in each duplex comes when the barcode
      # stays the same while both the order and mate switch). Process the duplex and start
      # a new one. If the barcode is the same, we're in the same duplex, but we've switched strands.
      if this_barcode != barcode or not (this_order != order and this_mate != mate):
        # sys.stderr.write('New duplex:  {}, {}, {}\n'.format(this_barcode, this_order, this_mate))
        process_duplex(duplex, barcode, workers=workers, stats=stats, **static)
        duplex = collections.OrderedDict()
      # else:
      #   sys.stderr.write('Same duplex: {}, {}, {}\n'.format(this_barcode, this_order, this_mate))
      barcode = this_barcode
      order = this_order
      mate = this_mate
      family = []
    read = {'name': name, 'seq':seq, 'qual':qual}
    family.append(read)
    all_reads += 1
  # Process the last family.
  duplex[(order, mate)] = family
  process_duplex(duplex, barcode, workers=workers, stats=stats, **static)

  if args.processes > 1:
    close_workers(workers)
    compile_results(workers)
    delete_tempfiles(workers)

  if args.sscs_file:
    static['sscs_fh'].close()
  if infile is not sys.stdin:
    infile.close()

  if not args.stats_file:
    return

  # Final stats on the run.
  logging.info('Processed {} reads and {} duplexes.'
               .format(all_reads, stats['runs']))
  per_read = stats['time'] / stats['reads']
  per_run = stats['time'] / stats['runs']
  logging.info('{:0.3f}s per read, {:0.3f}s per run.'.format(per_read, per_run))


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


def gather_args(args, infile, excluded_flags={'-S', '--slurm'},
                excluded_args={'-p', '--processes', '-l', '--log', '-s', '--sscs-file'}):
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


def delegate(worker, duplex, barcode):
  """Send a family to a worker process."""
  for (order, mate), family in duplex.items():
    for read in family:
      line = '{}\t{}\t{}\t{name}\t{seq}\t{qual}\n'.format(barcode, order, mate, **read)
    if family:
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


def process_duplex(duplex, barcode, workers=None, stats=None, incl_sscs=False, sscs_fh=None,
                   processes=1, min_reads=1, qual_thres=' '):
  stats['families'] += 1
  # Are we the controller process or a worker?
  if processes > 1:
    i = stats['families'] % len(workers)
    worker = workers[i]
    delegate(worker, duplex, barcode)
    return
  # We're a worker. Actually process the family.
  start = time.time()
  consensi = []
  reads_per_strand = []
  duplex_mate = None
  for (order, mate), family in duplex.items():
    reads = len(family)
    if reads < min_reads:
      continue
    # The mate number for the duplex consensus. It's arbitrary, but all that matters is that the
    # two mates have different numbers. This system ensures that:
    # Mate 1 is from the consensus of ab/1 and ba/2 families, while mate 2 is from ba/1 and ab/2.
    if (order == 'ab' and mate == 1) or (order == 'ba' and mate == 2):
      duplex_mate = 1
    else:
      duplex_mate = 2
    seqs = [read['seq'] for read in family]
    quals = [read['qual'] for read in family]
    consensi.append(consensus.get_consensus(seqs, quals, qual_thres=qual_thres))
    reads_per_strand.append(reads)
  assert len(consensi) <= 2
  if sscs_fh:
    for cons, (order, mate), reads in zip(consensi, duplex.keys(), reads_per_strand):
      sscs_fh.write('>{bar}.{order}.{mate} {reads}\n'.format(bar=barcode, order=order, mate=mate,
                                                             reads=reads))
      sscs_fh.write(cons+'\n')
  if len(consensi) == 1 and incl_sscs:
    print_duplex(consensi[0], barcode, duplex_mate, reads_per_strand)
  elif len(consensi) == 2:
    align = swalign.smith_waterman(*consensi)
    #TODO: log error & return if len(align.target) != len(align.query)
    cons = consensus.build_consensus_duplex_simple(align.target, align.query)
    print_duplex(cons, barcode, duplex_mate, reads_per_strand)
  elapsed = time.time() - start
  logging.info('{} sec for {} reads.'.format(elapsed, sum(reads_per_strand)))
  if stats and len(consensi) > 0:
    stats['time'] += elapsed
    stats['reads'] += sum(reads_per_strand)
    stats['runs'] += 1


def print_duplex(cons, barcode, mate, reads_per_strand, outfile=sys.stdout):
  header = '>{bar}.{mate} {reads}'.format(bar=barcode, mate=mate,
                                          reads='/'.join(map(str, reads_per_strand)))
  outfile.write(header+'\n')
  outfile.write(cons+'\n')


def read_fasta(fasta, is_file=True):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences. Each is a dict of 'name' and 'seq'.
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  seq_lines = []
  seq_name = None
  if is_file:
    with open(fasta) as fasta_file:
      fasta_lines = fasta_file.readlines()
  else:
    fasta_lines = fasta.splitlines()
  for line in fasta_lines:
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
