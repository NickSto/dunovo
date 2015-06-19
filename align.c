#include <stdio.h>
#include <string.h>
#include <ctype.h>

int _test_match(char *seq1, int start1, char *seq2, int start2);

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

#define NAIVE_TEST_WINDOW 6
#define NAIVE_TEST_THRES 0.80
#define NAIVE_TEST_MIN 2
#define NAIVE_WINDOW 10
#define NAIVE_THRES 0.80

// A naive algorithm for aligning two sequences which are expected to be very similar to each other
// and already nearly aligned.
void naive2(char *seq1, char *seq2) {
  int i = 0;
  int j = 0;
  int matches = 0;
  while (seq1[i] != 0 && seq2[j] != 0) {
    // Match?
    printf("%c %c | i %d j %d\n", seq1[i], seq2[j], i, j);
    if (seq1[i] == seq2[j]) {
      matches++;
      i++;
      j++;
      continue;
    }
    printf("mismatch!\n");
    // Mismatch. Start adding gaps until the mismatches go away.
    int new_i = i;
    int new_j = j;
    int gap_seq = 0;
    int success;
    while (1) {
      if (seq1[new_i] == 0 && seq2[new_j] == 0) {
        break;
      }
      success = _test_match(seq1, new_i, seq2, j);
      if (success) {
        gap_seq = 1;
        break;
      }
      if (seq1[new_i] != 0) {
        new_i++;
      }
      success = _test_match(seq1, i, seq2, new_j);
      if (success) {
        gap_seq = 2;
        break;
      }
      if (seq2[new_j] != 0) {
        new_j++;
      }
    }
    // Which sequence are we putting the gap in?
    if (gap_seq == 0) {
      printf("No good gap found. new_i: %d, new_j: %d\n", new_i, new_j);
      // No good gap found.
    } else if (i == new_i && j == new_j) {
      printf("No gap required.\n");
    } else if (gap_seq == 1) {
      printf("%dbp Gap in seq1 at %d.\n", new_i-i, i);
      i = new_i;
    } else if (gap_seq == 2) {
      printf("%dbp Gap in seq1 at %d.\n", new_j-j, j);
      j = new_j;
    }
  i++;
  j++;
  }
}

// Check if the few bases starting at start1 and start2 in seq1 and seq2, respectively, align with
// few mismatches. The number of bases checked is NAIVE_TEST_WINDOW, and they must have a match
// percentage greater than NAIVE_TEST_THRES. Also, the amount of sequence left to compare must be
// more than NAIVE_TEST_MIN.
int _test_match(char *seq1, int start1, char *seq2, int start2) {
  int matches = 0;
  int total = 0;
  char base1, base2;
  int i;
  for (i = 0; i < NAIVE_TEST_WINDOW-1; i++) {
    base1 = seq1[start1+i];
    base2 = seq2[start2+i];
    if (base1 == 0 || base2 == 0) {
      break;
    }
    if (base1 == base2) {
      matches++;
    }
    total++;
  }
  return total > NAIVE_TEST_MIN && (double)matches/total > NAIVE_TEST_THRES;
}
