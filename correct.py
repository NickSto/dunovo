#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import argparse

ARG_DEFAULTS = {'sam':sys.stdin, 'qual':20, 'pos':2, 'dist':1}
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

  args = parser.parse_args(argv[1:])

  # Group reads (barcodes) into sets of reads linked by alignments to each other.
  # Each group of reads is a set() of read names. "groups" is a list of these sets.
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
    # If the barcode is aligned to itself, we don't care about that.
    if read_name == ref_name:
      continue
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
          sys.stderr.write('Invalid NM tag "{}" on line {}\n.'.format(tag, line_num))
          raise
        break
    assert nm is not None, line_num
    if nm > args.dist:
      continue
    # It's a good alignment. Add it to a group.
    group_i_ref = member_to_group.get(ref_name)
    group_i_read = member_to_group.get(read_name)
    if group_i_ref is None:
      group_i = group_i_read
    elif group_i_read is None:
      group_i = group_i_ref
    elif group_i_ref == group_i_read:
      group_i = group_i_ref
    else:
      # Both the read and ref are already in a group, but they're different groups. We need to join
      # them.
      join_groups(groups, member_to_group, ref_name, read_name)
      group_i = member_to_group[ref_name]
    if group_i is None:
      new_group = set((ref_name, read_name))
      groups.append(new_group)
      group_i = len(groups) - 1
      member_to_group[ref_name] = group_i
      member_to_group[read_name] = group_i
    else:
      group = groups[group_i]
      group.add(ref_name)
      group.add(read_name)


def join_groups(groups, member_to_group, ref_name, read_name):
  group_i_ref = member_to_group[ref_name]
  group_i_read = member_to_group[read_name]
  group_ref = groups[group_i_ref]
  group_read = groups[group_i_read]
  # Get a union of all elements in both groups.
  group_union = group_ref.union(group_read)
  # Put the new, union group back in place of the ref group, and make that the canonical group.
  groups[group_i_ref] = group_union
  groups[group_i_read] = None
  member_to_group[ref_name] = group_i_ref
  member_to_group[read_name] = group_i_ref


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
