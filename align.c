#include <stdio.h>
#include <string.h>

// Take an existing alignment and consensus and compute the number of differences between each
// sequence and the consensus.
int *get_diffs(int n_seqs, char *seqs[], char *cons) {
  int diffs[n_seqs];
  int i, j;
  char *seq;
  // Loop through the sequences in the alignment.
  for (i = 0; i < n_seqs; i++) {
    j = 0;
    seq = seqs[i];
    diffs[i] = 0;
    // Compare each base of the sequence to the consensus.
    while (seq[j] != 0 && cons[j] != 0) {
      if (seq[j] != cons[j]) {
        diffs[i]++;
      }
      j++;
    }
  }
  return diffs;
}

