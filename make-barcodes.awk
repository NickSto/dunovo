# The awk code that transforms the one-line fastq record pair into the output that can be sorted
# by barcode.
# The input must have 8 columns: the 4 FASTQ lines for both reads in a read pair.
# Output columns:
#   1: the barcode, put into a canonical form
#   2: the order of the barcode halves ("ab" or "ba")
#   3: read1 name
#   4: sequence of read 1, minus the 12bp barcode and 5bp constant sequence
#   5: read1 quality scores, minus the same first 17bp
#   6: read2 name
#   7: sequence of read 2, minus the first 17bp
#   8: read2 quality scores, minus the first 17bp
# The canonical form of the barcode is composed of the first 12bp of each read. The 12bp from the
# first read is the "alpha" and the one from the second is the "beta". The barcode is formed by
# concatenating them in an order determined by summing the ASCII values of each 12bp sequence, then
# putting the one with a lower sum first. If the alpha is first, the value of the second column (the
# order) is "ab", and "ba" if the beta is first.

BEGIN {
  FS = "\t"; OFS = "\t";
  if (! BAR_LEN) {
    BAR_LEN = 24;
  }
  if (! INVARIANT) {
    INVARIANT = 5;
  }
  HALF = int(BAR_LEN/2);
  # Fill an array mapping each character to its ASCII value.
  for (i = 0; i < 128; i++) {
    char = sprintf("%c", i);
    ORD[char] = i;
  }
}

function get_sum(seq) {
  sum = 0;
  for (i = 1; i <= length(seq); i++) {
    sum += ORD[substr(seq, i, 1)];
  }
  return sum;
}

{
  alpha = substr($2, 1, HALF);
  beta = substr($6, 1, HALF);
  alpha_sum = get_sum(alpha);
  beta_sum = get_sum(beta);
  if (alpha_sum < beta_sum) {
    barcode = alpha beta;
    order = "ab";
  } else {
    barcode = beta alpha;
    order = "ba";
  }
  name1 = $1;
  name2 = $5;
  seq1 = substr($2, HALF + INVARIANT + 1);
  seq2 = substr($6, HALF + INVARIANT + 1);
  qual1 = substr($4, HALF + INVARIANT + 1);
  qual2 = substr($8, HALF + INVARIANT + 1);
  print barcode, order, name1, seq1, qual1, name2, seq2, qual2;
}
