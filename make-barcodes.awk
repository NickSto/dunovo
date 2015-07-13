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
  read1 = substr($2, HALF + INVARIANT + 1);
  read2 = substr($6, HALF + INVARIANT + 1);
  # $1 = read1 name | $4 = read1 quality | $5 = read2 name | $8 = read2 quality
  print barcode, order, $1, read1, $4, $5, read2, $8;
}

