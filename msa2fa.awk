# A quick script to convert the .msa.tsv output of sscs.py back into FASTA format.

BEGIN {
  FS = "\t";
  OFS = "\t";
}

$2 == "CONSENSUS" {
  if ($1 == last) {
    mate = 2;
  } else {
    mate = 1;
  }
  printf(">%s.%d:%d\n", $1, mate, pairs);
  print $3;
  pairs = 0;
  last = $1;
}

$2 != "CONSENSUS" {
  pairs++;
}