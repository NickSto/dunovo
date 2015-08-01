import os
import ctypes

script_dir = os.path.dirname(os.path.realpath(__file__))
consensus = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'consensusc.so'))
consensus.get_consensus.restype = ctypes.c_char_p


# N.B.: The quality scores must be aligned with their accompanying sequences.
def get_consensus(align, quals=[], cons_thres=-1.0, qual_thres=' ', gapped=False):
  cons_thres_c = ctypes.c_double(cons_thres)
  qual_thres_c = ctypes.c_char(qual_thres)
  n_seqs = len(align)
  if gapped:
    gapped_c = 1
  else:
    gapped_c = 0
  assert not quals or len(quals) == n_seqs, 'Different number of sequences and quals.'
  seq_len = None
  for seq in (align + quals):
    if seq_len is None:
      seq_len = len(seq)
    else:
      assert seq_len == len(seq), 'All sequences in the alignment must be the same length.'
  align_c = (ctypes.c_char_p * n_seqs)()
  for i, seq in enumerate(align):
    align_c[i] = ctypes.c_char_p(seq)
  quals_c = (ctypes.c_char_p * n_seqs)()
  for i, qual in enumerate(quals):
    quals_c[i] = ctypes.c_char_p(qual)
  if not quals:
    quals_c = 0
  return consensus.get_consensus(align_c, quals_c, n_seqs, seq_len, cons_thres_c, qual_thres_c,
                                 gapped_c)
