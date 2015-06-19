#include <stdio.h>
#include <string.h>
#include <ctype.h>

// Take an existing alignment and consensus and compute the number of differences between each
// sequence and the consensus.
int *get_diffs(char *cons, char *seqs[], int n_seqs) {
  int diffs[n_seqs];
  int i = 0;
  // Uppercase the consensus.
  while (cons[i] != 0) {
    cons[i] = toupper(cons[i]);
    i++;
  }
  int j;
  char *seq;
  // Loop through the sequences in the alignment.
  for (i = 0; i < n_seqs; i++) {
    j = 0;
    seq = seqs[i];
    diffs[i] = 0;
    // Compare each base of the sequence to the consensus.
    while (seq[j] != 0 && cons[j] != 0) {
      if (toupper(seq[j]) != cons[j]) {
        diffs[i]++;
      }
      j++;
    }
  }
  return diffs;
}
