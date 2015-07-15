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
  print ">" $1 "." mate;
  print $3;
  last = $1;
}
