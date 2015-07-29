import os
import ctypes

script_dir = os.path.dirname(os.path.realpath(__file__))
consensus = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'consensusc.so'))
consensus.get_consensus.restype = ctypes.c_char_p


def get_consensus(alignment, cons_thres=-1.0, qualities=None, qual_thres=' '):
  cons_thres_c = ctypes.c_double(cons_thres)
  qual_thres_c = ctypes.c_char(qual_thres)
  n_seqs = len(alignment)
  assert qualities is None or len(qualities) == n_seqs, 'Different number of sequences and qualities.'
  seq_len = None
  for seq in alignment:
    if seq_len is None:
      seq_len = len(seq)
    else:
      assert seq_len == len(seq), 'All sequences in the alignment must be the same length.'
  align_c = (ctypes.c_char_p * n_seqs)()
  for i, seq in enumerate(alignment):
    align_c[i] = ctypes.c_char_p(seq)
  if qualities:
    qual_c = (ctypes.c_char_p * n_seqs)()
    for i, qual in enumerate(qualities):
      qual_c[i] = ctypes.c_char_p(qual)
  else:
    qual_c = 0
  return consensus.get_consensus(align_c, qual_c, n_seqs, seq_len, cons_thres_c, qual_thres_c)
