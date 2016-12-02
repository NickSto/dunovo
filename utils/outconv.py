#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import argparse

ARG_DEFAULTS = {'input':sys.stdin}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='consensuses.fa', type=argparse.FileType('r'),
    help='Interleaved consensus reads (DCS or SSCS) from dunovo.py.')
  parser.add_argument('-1', '--out1', metavar='out_1.fa', type=argparse.FileType('w'), required=True,
    help='Output filename for mate 1 reads. CAUTION: Will overwrite any existing file.')
  parser.add_argument('-2', '--out2', metavar='out_2.fa', type=argparse.FileType('w'), required=True,
    help='Output filename for mate 2 reads. CAUTION: Will overwrite any existing file.')

  args = parser.parse_args(argv[1:])

  """
SSCS:
>AAAAAAAAAAAACTAAAATACAAA.ab.1 15
ACTGATGAAAAGGCTGTTATTGTATCTGATGTGTAGTGTATGGCTAAGAAAAGACCTGTAATGATTTGGACTATTAGGCAGACTCCTAGAAGGGACCCAA
>AAAAAAAAAAAACTAAAATACAAA.ab.2 15
ATCCAAACACAACCAACATCCCCCCTAAATAAATTAAAAAAACTATTAAACCTAAAAACGATCCACCAAACCCTAAAACCATTAAACAACCAACAAACCC
>AAAAAAAAAAAAGACCACGTTTCT.ab.1 15
TGGAATTCAGCCTACTAGCAATTATCCCCATACTAATCAACAAAAAAAACCCACGATCAACTGAAGCAGCAACAAAATACTTCGTCACACAAGCAACAGC
>AAAAAAAAAAAAGACCACGTTTCT.ba.2 18
TGGAATTCAGCCTACTAGCAATTATCCCCATACTAATCAACAAAAAAAACCCACGATCAACTGAAGCAGCAACAAAATACTTCGTCACACAAGCAACAGC
>AAAAAAAAAAAAGACCACGTTTCT.ab.2 15
GTAGAGTTGAGTAGCGGGTAAATTTGAATTAAAATTGATAGGGGAGCAATTTTTTGTCATGTAAGAAGAATAAGTCCTATGTGCAGTGGGATCCCTTGAG
>AAAAAAAAAAAAGACCACGTTTCT.ba.1 18
GTAGAGTTGAGTAGCGGGTAAATTTGAATTAAAATTGATAGGGGAGCAATTTTTTGTCATGTAAGAAGAATAAGTCCTATGTGCAGTGGGATCCCTTGAG

DCS:
>AAAAAAAACACCAAATACGCCTAC.1 28-39
GCTATGAATATAGGGGCTGTAAGAATAATATAGATTATGAGGTTGAGTAGAGTGAGGGATGGGTTGTAAGGAAGAATTGCTAATATTCATCCTATGTGGG
>AAAAAAAACACCAAATACGCCTAC.2 28-39
TACTTCGTCACACAAGCAACAGCCTCAATAATTATCCTCCTGGCCATCGTACTCAACTATAAACAACTAGGAACATGAATATTTCAACAACAAACAAACG
  """

  intype = None
  line_num = 0
  sscs_buffer = {}
  for line_raw in args.input:
    line_num += 1
    line = line_raw.rstrip('\r\n')
    if line.startswith('>'):
      # We're in a header line.
      barcode, strand, mate, famsizes = parse_header(line)
      if intype is None:
        if strand is None:
          intype = 'DCS'
        else:
          intype = 'SSCS'
    else:
      # We're in a sequence line.
      seq = line
      if intype == 'DCS':
        if mate == '1':
          args.out1.write('>{} {}\n{}\n'.format(barcode, famsizes, seq))
        elif mate == '2':
          args.out2.write('>{} {}\n{}\n'.format(barcode, famsizes, seq))
        else:
          fail('Error: Invalid mate "{}" on line {}.'.format(mate, line_num))
      elif intype == 'SSCS':
        # Need to get SSCS in correct order.
        # Collect reads until they're properly paired, then print both.
        if (barcode, strand) in sscs_buffer:
          stored_read = sscs_buffer[(barcode, strand)]
          new_read = {'mate':mate, 'famsizes':famsizes, 'seq':seq}
          if stored_read['mate'] == '1' and mate == '2':
            read1, read2 = stored_read, new_read
          elif stored_read['mate'] == '2' and mate == '1':
            read1, read2 = new_read, stored_read
          args.out1.write('>{}.{} {famsizes}\n{seq}\n'.format(barcode, strand, **read1))
          args.out2.write('>{}.{} {famsizes}\n{seq}\n'.format(barcode, strand, **read2))
          del sscs_buffer[(barcode, strand)]
        else:
          sscs_buffer[(barcode, strand)] = {'mate':mate, 'famsizes':famsizes, 'seq':seq}


def parse_header(header):
  read_id, famsizes = header.split()
  id_fields = read_id[1:].split('.')
  barcode = id_fields[0]
  strand = None
  mate = id_fields[-1]
  if len(id_fields) == 3:
    strand = id_fields[1]
  return barcode, strand, mate, famsizes


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
