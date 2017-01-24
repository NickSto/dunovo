#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import gzip
import logging
import argparse
import collections

VERBOSE = (logging.DEBUG+logging.INFO)//2
ARG_DEFAULTS = {'sam':sys.stdin, 'qual':20, 'pos':2, 'dist':1, 'log':sys.stderr,
                'volume':logging.WARNING}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Correct barcodes using an alignment of all barcodes to themselves. Reads the
alignment in SAM format and corrects the barcodes in an input "families" file (the output of
make-barcodes.awk). It will print the "families" file to stdout with barcodes (and orders)
corrected."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('families', type=open_as_text_or_gzip,
    help='The sorted output of make-barcodes.awk. The important part is that it\'s a tab-delimited '
         'file with at least 2 columns: the barcode sequence and order, and it must be sorted in '
         'the same order as the "reads" in the SAM file.')
  parser.add_argument('reads', type=open_as_text_or_gzip,
    help='The fasta/q file given to the aligner. Used to get barcode sequences from read names.')
  parser.add_argument('sam', type=argparse.FileType('r'), nargs='?',
    help='Barcode alignment, in SAM format. Omit to read from stdin. The read names must be '
         'integers, representing the (1-based) order they appear in the families file.')
  parser.add_argument('-P', '--prepend', action='store_true',
    help='Prepend the corrected barcodes and orders to the original columns.')
  parser.add_argument('-d', '--dist', type=int,
    help='NM edit distance threshold. Default: %(default)s')
  parser.add_argument('-m', '--mapq', type=int,
    help='MAPQ threshold. Default: %(default)s')
  parser.add_argument('-p', '--pos', type=int,
    help='POS tolerance. Alignments will be ignored if abs(POS - 1) is greater than this value. '
         'Set to greater than the barcode length for no threshold. Default: %(default)s')
  parser.add_argument('-L', '--tag-len', type=int,
    help='Length of each half of the barcode. If not given, it will be determined from the first '
         'barcode in the families file.')
  parser.add_argument('-o', '--reorder', action='store_true')
  parser.add_argument('--limit', type=int,
    help='Limit the number of lines that will be read from each input file, for testing purposes.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-i', '--info', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=VERBOSE)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG,
    help='Print debug messages (very verbose).')

  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  logging.info('Reading the fasta/q to map read names to barcodes..')
  names_to_barcodes = map_names_to_barcodes(args.reads, args.limit)

  logging.info('Reading the SAM to build groups from SAM alignment..')
  correction_table = make_correction_table(args.sam, names_to_barcodes, args.pos, args.mapq,
                                           args.dist, args.limit)

  if args.reorder:
    logging.info('Reading the families.tsv to get the orders of the barcodes..')
    orders = get_orders(args.families, args.limit)
  else:
    orders = None

  logging.info('Reading the families.tsv again to print corrected output..')
  families = open_as_text_or_gzip(args.families.name)
  print_corrected_output(families, correction_table, orders, args.prepend, args.limit)


def detect_format(reads_file, max_lines=7):
  """Detect whether a file is a fastq or a fasta, based on its content."""
  fasta_votes = 0
  fastq_votes = 0
  line_num = 0
  for line in reads_file:
    line_num += 1
    if line_num % 4 == 1:
      if line.startswith('@'):
        fastq_votes += 1
      elif line.startswith('>'):
        fasta_votes += 1
    elif line_num % 4 == 3:
      if line.startswith('+'):
        fastq_votes += 1
      elif line.startswith('>'):
        fasta_votes += 1
    if line_num >= max_lines:
      break
  reads_file.seek(0)
  if fasta_votes > fastq_votes:
    return 'fasta'
  elif fastq_votes > fasta_votes:
    return 'fastq'
  else:
    return None


def read_fastaq(reads_file):
  filename = reads_file.name
  if filename.endswith('.fa') or filename.endswith('.fasta'):
    format = 'fasta'
  elif filename.endswith('.fq') or filename.endswith('.fastq'):
    format = 'fastq'
  else:
    format = detect_format(reads_file)
  if format == 'fasta':
    return read_fasta(reads_file)
  elif format == 'fastq':
    return read_fastq(reads_file)


def read_fasta(reads_file):
  """Read a FASTA file, yielding read names and sequences.
  NOTE: This assumes sequences are only one line!"""
  line_num = 0
  for line_raw in reads_file:
    line = line_raw.rstrip('\r\n')
    line_num += 1
    if line_num % 2 == 1:
      assert line.startswith('>'), line
      read_name = line[1:]
    elif line_num % 2 == 0:
      read_seq = line
      yield read_name, read_seq


def read_fastq(reads_file):
  """Read a FASTQ file, yielding read names and sequences.
  NOTE: This assumes sequences are only one line!"""
  line_num = 0
  for line in reads_file:
    line_num += 1
    if line_num % 4 == 1:
      assert line.startswith('@'), line
      read_name = line[1:].rstrip('\r\n')
    elif line_num % 4 == 2:
      read_seq = line.rstrip('\r\n')
      yield read_name, read_seq


def map_names_to_barcodes(reads_file, limit=None):
  """Map barcode names to their sequences."""
  names_to_barcodes = {}
  read_num = 0
  for read_name, read_seq in read_fastaq(reads_file):
    read_num += 1
    if limit is not None and read_num > limit:
      break
    try:
      name = int(read_name)
    except ValueError:
      logging.critical('non-int read name "{}"'.format(name))
      raise
    names_to_barcodes[name] = read_seq
  reads_file.close()
  return names_to_barcodes


def parse_alignment(sam_file, pos_thres, mapq_thres, dist_thres):
  """Parse the SAM file and yield reads that pass the filters.
  Returns (read_name, ref_name)."""
  line_num = 0
  for line in sam_file:
    line_num += 1
    if line.startswith('@'):
      logging.debug('Header line ({})'.format(line_num))
      continue
    fields = line.split('\t')
    logging.debug('read {} -> ref {} (read seq {}):'.format(fields[2], fields[0], fields[9]))
    try:
      read_name = int(fields[0])
      ref_name = int(fields[2])
    except ValueError:
      if fields[2] == '*':
        logging.debug('\tRead unmapped (reference == "*")')
        continue
      else:
        logging.error('Non-integer read name(s) on line {}: "{}", "{}".'
                      .format(line_num, read_name, ref_name))
        raise
    # Apply alignment quality filters.
    try:
      flags = int(fields[1])
      pos = int(fields[3])
      mapq = int(fields[4])
    except ValueError:
      logging.warn('\tNon-integer flag ({}), pos ({}), or mapq ({})'
                   .format(fields[1], fields[3], fields[4]))
      continue
    if flags & 4:
      logging.debug('\tRead unmapped (flag & 4 == True)')
      continue
    if abs(pos - 1) > pos_thres:
      logging.debug('\tAlignment failed pos filter: abs({} - 1) > {}'.format(pos, pos_thres))
      continue
    if mapq < mapq_thres:
      logging.debug('\tAlignment failed mapq filter: {} > {}'.format(mapq, mapq_thres))
      continue
    nm = None
    for tag in fields[11:]:
      if tag.startswith('NM:i:'):
        try:
          nm = int(tag[5:])
        except ValueError:
          logging.error('Invalid NM tag "{}" on line {}.'.format(tag, line_num))
          raise
        break
    assert nm is not None, line_num
    if nm > dist_thres:
      logging.debug('\tAlignment failed NM distance filter: {} > {}'.format(nm, dist_thres))
      continue
    yield read_name, ref_name
  sam_file.close()


def make_correction_table(sam_file, names_to_barcodes, pos_thres, mapq_thres, dist_thres, limit=None):
  """Make a table mapping original barcode numbers to correct barcode numbers."""
  correction_table = {}
  # Maps correct barcode numbers to sets of original barcodes (includes correct ones).
  reverse_table = collections.defaultdict(set)
  line_num = 0
  for read_name, ref_name in parse_alignment(sam_file, pos_thres, mapq_thres, dist_thres):
    line_num += 1
    if limit is not None and line_num > limit:
      break
    # Skip self-alignments.
    if ref_name == read_name:
      continue
    correct = names_to_barcodes[ref_name]
    original = names_to_barcodes[read_name]
    # Check if the "correct" barcode has already been corrected.
    if correct in correction_table:
      corrected_correct = correction_table[correct]
      correction_table[original] = corrected_correct
      group = reverse_table[corrected_correct]
      group.add(correct)
      group.add(original)
      group.add(corrected_correct)  # Not necessary, but just in case.
      # And that's all we need to do.
      continue
    # Check if this barcode has already been corrected to something else.
    #TODO: What if the "correct" barcode has already been corrected
    #      AND the original has already been corrected to something else?
    if original in correction_table:
      # If so, override that correction.
      # First, change the mapping in reverse_table.
      old_correct = correction_table[original]
      group = reverse_table[old_correct]
      reverse_table[correct] = group
      logging.debug('Read {} already corrected to {}. Changing correction to {}.'
                    .format(original, old_correct, correct))
      # Then, change the mapping for each member of the group.
      for member in group:
        correction_table[member] = correct
    # Create this mapping.
    correction_table[original] = correct
    group = reverse_table[correct]
    group.add(original)
    group.add(correct)
  return correction_table


def get_orders(families_file, limit=None):
  """Map barcode sequences to their order strings ("ab" or "ba")."""
  orders = {}
  line_num = 0
  for line in families_file:
    line_num += 1
    if limit is not None and line_num > limit:
      break
    fields = line.rstrip('\r\n').split('\t')
    barcode = fields[0]
    order = fields[1]
    orders[barcode] = order
  families_file.close()
  return orders


def print_corrected_output(families_file, corrections, orders, prepend=False, tag_len=None, limit=None):
  # Determine barcode tag length if not given.
  if tag_len is None:
    tag_len = len(corrections.keys()[1])//2
  line_num = 0
  barcode_num = 0
  barcode_last = None
  correct_order = None
  corrected = {'reads':0, 'barcodes':0, 'orders':0}
  reads = [0, 0]
  corrections_in_this_family = 0
  for line in families_file:
    line_num += 1
    if limit is not None and line_num > limit:
      break
    fields = line.rstrip('\r\n').split('\t')
    raw_barcode = fields[0]
    raw_order = fields[1]
    if raw_barcode != barcode_last:
      # We just started a new family.
      barcode_num += 1
      family_info = '{}\t{}\t{}'.format(barcode_last, reads[0], reads[1])
      if corrections_in_this_family:
        corrected['reads'] += corrections_in_this_family
        corrected['barcodes'] += 1
        if orders and raw_order != correct_order:
          corrected['orders'] += 1
        family_info += '\tCORRECTED!'
      else:
        family_info += '\tuncorrected'
      logging.log(VERBOSE, family_info)
      reads = [0, 0]
      corrections_in_this_family = 0
      barcode_last = raw_barcode
    if raw_order == 'ab':
      reads[0] += 1
    elif raw_order == 'ba':
      reads[1] += 1
    if raw_barcode in corrections:
      correct_barcode = corrections[raw_barcode]
      if orders:
        correct_order = orders[correct_barcode]
      corrections_in_this_family += 1
    else:
      correct_barcode = raw_barcode
      correct_order = raw_order
    if prepend:
      fields.insert(0, correct_barcode)
      fields.insert(0, correct_order)
    else:
      fields[0] = correct_barcode
      fields[1] = correct_order
    print(*fields, sep='\t')
  families_file.close()
  if corrections_in_this_family:
    corrected['reads'] += corrections_in_this_family
    corrected['barcodes'] += 1
  logging.info('Corrected {barcodes} barcodes on {reads} read pairs, with {orders} order reversals.'
               .format(**corrected))


def open_as_text_or_gzip(path):
  """Return an open file-like object reading the path as a text file or a gzip file, depending on
  which it looks like."""
  if detect_gzip(path):
    return gzip.open(path, 'r')
  else:
    return open(path, 'rU')


def detect_gzip(path):
  """Return True if the file looks like a gzip file: ends with .gz or contains non-ASCII bytes."""
  ext = os.path.splitext(path)[1]
  if ext == '.gz':
    return True
  elif ext in ('.txt', '.tsv', '.csv'):
    return False
  with open(path) as fh:
    is_not_ascii = detect_non_ascii(fh.read(100))
  if is_not_ascii:
    return True


def detect_non_ascii(bytes, max_test=100):
  """Return True if any of the first "max_test" bytes are non-ASCII (the high bit set to 1).
  Return False otherwise."""
  for i, char in enumerate(bytes):
    # Is the high bit a 1?
    if ord(char) & 128:
      return True
    if i >= max_test:
      return False
  return False


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


if __name__ == '__main__':
  sys.exit(main(sys.argv))
