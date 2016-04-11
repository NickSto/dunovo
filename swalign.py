import os
import ctypes

script_dir = os.path.dirname(os.path.realpath(__file__))
swalign = ctypes.cdll.LoadLibrary(os.path.join(script_dir, 'libswalign.so'))


# C struct for ctypes
class SeqPairC(ctypes.Structure):
  _fields_ = [
    ('a', ctypes.c_char_p),
    ('alen', ctypes.c_uint),
    ('b', ctypes.c_char_p),
    ('blen', ctypes.c_uint),
  ]


# C struct for ctypes
class AlignC(ctypes.Structure):
  _fields_ = [
    ('seqs', ctypes.POINTER(SeqPairC)),
    ('start_a', ctypes.c_int),
    ('start_b', ctypes.c_int),
    ('end_a', ctypes.c_int),
    ('end_b', ctypes.c_int),
    ('matches', ctypes.c_int),
    ('score', ctypes.c_double),
  ]


# The Python version
class Align(object):
  def __init__(self, align_c):
    self.target = align_c.seqs.contents.a
    self.query = align_c.seqs.contents.b
    self.start_target = align_c.start_a
    self.start_query = align_c.start_b
    self.end_target = align_c.end_a
    self.end_query = align_c.end_b
    self.matches = align_c.matches
    self.score = align_c.score


# Initialize functions (define types).
swalign.smith_waterman.restype = ctypes.POINTER(AlignC)
swalign.revcomp.restype = ctypes.c_char_p


def smith_waterman(target, query):
  seq_pair = SeqPairC(target, len(target), query, len(query))
  align_c = swalign.smith_waterman(ctypes.pointer(seq_pair), 1).contents
  return Align(align_c)


def revcomp(seq):
  """WARNING: This will alter the input string in-place!"""
  swalign.revcomp(seq)
