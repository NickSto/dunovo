#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import logging
import argparse

ARG_DEFAULTS = {'sam':sys.stdin, 'qual':20, 'pos':2, 'dist':1, 'log':sys.stderr,
                'log_level':logging.ERROR}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('sam', type=argparse.FileType('r'), nargs='?',
    help='Barcode alignment, in SAM format. Omit to read from stdin.')
  parser.add_argument('-d', '--dist', type=int,
    help='NM edit distance threshold. Default: %(default)s')
  parser.add_argument('-q', '--qual', type=int,
    help='MAPQ threshold. Default: %(default)s')
  parser.add_argument('-p', '--pos', type=int,
    help='POS tolerance. Alignments will be ignored if abs(POS - 1) is greater than this value. '
         'Set to greater than the barcode length for no threshold. Default: %(default)s')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-D', '--debug', dest='log_level', action='store_const', const=logging.DEBUG,
    help='Print debug messages (very verbose).')

  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.log_level, format='%(message)s')
  tone_down_logger()

  # Group reads (barcodes) into sets of reads linked by alignments to each other.
  # Each group of reads is a dict mapping read names to read sequences. "groups" is a list of these
  # dicts.
  groups = []
  # "member_to_group" is a dict mapping read names to their group's index in "groups".
  member_to_group = {}
  line_num = 0
  for line in args.sam:
    line_num += 1
    if line.startswith('@'):
      continue
    fields = line.split('\t')
    read_name = fields[0]
    ref_name = fields[2]
    read_seq = fields[9]
    # Apply alignment quality filters.
    try:
      flags = int(fields[1])
      pos = int(fields[3])
      mapq = int(fields[4])
    except ValueError:
      continue
    if flags & 4:
      # Read unmapped.
      continue
    if abs(pos - 1) > args.pos:
      continue
    if mapq < args.qual:
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
    if nm > args.dist:
      continue
    logging.debug('Ref {}, read {}:'.format(ref_name, read_name))
    # It's a good alignment. Add it to a group.
    group_i_ref = member_to_group.get(ref_name)
    group_i_read = member_to_group.get(read_name)
    if group_i_ref is None:
      group_i = group_i_read
      logging.debug('\tgroup_i_ref  is None, using {}.'.format(group_i_read))
    elif group_i_read is None:
      group_i = group_i_ref
      logging.debug('\tgroup_i_read is None, using {}.'.format(group_i_ref))
    elif group_i_ref == group_i_read:
      group_i = group_i_ref
      logging.debug('\tgroup_i_ref == group_i_read ({}).'.format(group_i_ref))
    else:
      # Both the read and ref are already in a group, but they're different groups. We need to join
      # them.
      logging.debug('\tgroup_i_ref  is {}, group_i_read is {}.'.format(group_i_ref, group_i_read))
      logging.debug('\tJoining groups {} and {}.'.format(group_i_ref, group_i_read))
      join_groups(groups, member_to_group, ref_name, read_name)
      group_i = member_to_group[ref_name]
    if group_i is None:
      # No group exists yet. Create one and add it.
      new_group = {ref_name:None, read_name:read_seq}
      groups.append(new_group)
      group_i = len(groups) - 1
      member_to_group[ref_name] = group_i
      member_to_group[read_name] = group_i
      logging.debug('\tAdding new group ({}): {} and {}.'.format(group_i, ref_name, read_name))
    else:
      # Add these reads to the group.
      logging.debug('\tAdding {} and {} to group {}.'.format(ref_name, read_name, group_i))
      group = groups[group_i]
      if not group.get(ref_name):
        logging.debug('\t\t{} not in {}; adding it.'.format(ref_name, group_i))
        group[ref_name] = None
        member_to_group[ref_name] = group_i
      else:
        logging.debug('\t\t{} already in {}.'.format(ref_name, group_i))
      if not group.get(read_name):
        logging.debug('\t\t{} not in {}; adding it.'.format(read_name, group_i))
        group[read_name] = read_seq
        member_to_group[read_name] = group_i
      else:
        logging.debug('\t\t{} already in {}.'.format(read_name, group_i))

  for i, group in enumerate(groups):
    print('{}:'.format(i))
    if group:
      for name, seq in group.items():
        print('\t{}:\t{}'.format(name, seq))


def join_groups(groups, member_to_group, ref_name, read_name):
  group_i_ref = member_to_group[ref_name]
  group_i_read = member_to_group[read_name]
  logging.debug('\t\tmember_to_group[{}] == {}'.format(ref_name, group_i_ref))
  logging.debug('\t\tmember_to_group[{}] == {}'.format(read_name, group_i_read))
  group_ref = groups[group_i_ref]
  group_read = groups[group_i_read]
  # Get a union of all elements in both groups.
  group_union = group_ref
  for name, seq in group_read.items():
    # Only update the read if its name isn't present or its sequence is None.
    if not group_union.get(name):
      group_union[name] = seq
  # Put the new, union group back in place of the ref group, and make that the canonical group.
  groups[group_i_ref] = group_union
  groups[group_i_read] = None
  for name in group_union:
    member_to_group[name] = group_i_ref


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


if __name__ == '__main__':
  sys.exit(main(sys.argv))
