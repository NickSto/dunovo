import os
import ctypes

script_dir = os.path.dirname(os.path.realpath(__file__))
seqtools = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'seqtoolsc.so'))
seqtools.get_revcomp.restype = ctypes.c_char_p


def get_revcomp(seq):
  return seqtools.get_revcomp(seq)


def get_diffs_frac_simple(consensus, family):
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    c_family[i] = ctypes.c_char_p(seq)
  seqtools.get_diffs_frac_simple.restype = ctypes.POINTER(ctypes.c_double * len(c_family))
  diffs = seqtools.get_diffs_frac_simple(c_consensus, c_family, len(c_family))
  return diffs.contents


def get_diffs_frac_binned(consensus, family, bins):
  seq_len = None
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    if seq_len:
      if seq_len != len(seq):
        return None
    else:
      seq_len = len(seq)
    c_family[i] = ctypes.c_char_p(seq)
  double_array_pointer = ctypes.POINTER(ctypes.c_double * bins)
  seqtools.get_diffs_frac_binned.restype = ctypes.POINTER(double_array_pointer * len(c_family))
  diffs_binned_c = seqtools.get_diffs_frac_binned(c_consensus, c_family, len(c_family), seq_len, bins)
  diffs_binned = []
  for diffs_c in diffs_binned_c.contents:
    diffs_binned.append(diffs_c.contents)
  return diffs_binned
