#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct Gap {
  int seq;
  int coord;
  int length;
  struct Gap *next;
} Gap;

typedef struct Gaps {
  int length;
  struct Gap *root;
  struct Gap *tip;
} Gaps;

int _test_match(char *seq1, int start1, char *seq2, int start2);
void add_gap(Gaps *gaps, int seq, int coord, int length);
Gaps *make_gaps();
char *insert_gaps(Gaps *gaps, char *seq, int seq_num);

// Take an existing alignment and consensus and compute the number of differences between each
// sequence and the consensus.
int *get_diffs_simple(char *cons, char *seqs[], int n_seqs) {
  int *diffs = malloc(sizeof(int) * n_seqs);
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
  Gaps *gaps = make_gaps();
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
        gap_seq = 2;
        break;
      }
      if (seq1[new_i] != 0) {
        new_i++;
      }
      success = _test_match(seq1, i, seq2, new_j);
      if (success) {
        gap_seq = 1;
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
      printf("%dbp gap in seq1 at base %d.\n", new_j-j, j);
      add_gap(gaps, 1, j, new_j-j);
      j = new_j;
    } else if (gap_seq == 2) {
      printf("%dbp gap in seq2 at base %d.\n", new_i-i, i);
      add_gap(gaps, 2, i, new_i-i);
      i = new_i;
    }
    i++;
    j++;
  }

  char *new_seq1 = insert_gaps(gaps, seq1, 1);
  char *new_seq2 = insert_gaps(gaps, seq2, 2);
  printf("alignment:\n%s\n%s\n", new_seq1, new_seq2);
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

Gaps *make_gaps() {
  Gaps *gaps = malloc(sizeof(Gaps));
  gaps->root = 0;
  gaps->tip = 0;
  gaps->length = 0;
  return gaps;
}

void add_gap(Gaps *gaps, int seq, int coord, int length) {
  Gap *gap = malloc(sizeof(Gap));
  gap->next = 0;
  gap->seq = seq;
  gap->coord = coord;
  gap->length = length;
  if (gaps->root == 0) {
    gaps->root = gap;
  } else {
    gaps->tip->next = gap;
  }
  gaps->tip = gap;
  gaps->length++;
}

// Take gap information from the aligner and put them into the sequence string as "-" characters.
char *insert_gaps(Gaps *gaps, char *seq, int seq_num) {
  if (gaps->root == 0) {
    return seq;
  }

  // How long should the new sequence be?
  int extra_len = 0;
  Gap *gap = gaps->root;
  while (gap) {
    if (gap->seq == seq_num) {
      extra_len += gap->length;
    }
    gap = gap->next;
  }

  //TODO: Handle a situation with no gaps.
  int new_len = extra_len + strlen(seq) + 1;
  char *new_seq = malloc(sizeof(char) * new_len);
  int i = 0;
  int j = 0;
  gap = gaps->root;
  while (gap) {
    // Check that it's a gap in our sequence.
    if (gap->seq != seq_num) {
      gap = gap->next;
      continue;
    }
    // Copy verbatim all the sequence until the gap.
    while (i <= gap->coord) {
      new_seq[j] = seq[i];
      i++;
      j++;
    }
    // Add -'s the whole length of the gap.
    while (j < gap->coord + gap->length + 1) {
      new_seq[j] = '-';
      j++;
    }
    gap = gap->next;
  }
  // Fill in the end sequence.
  while (seq[i]) {
    new_seq[j] = seq[i];
    i++;
    j++;
  }
  new_seq[new_len-1] = 0;
  return new_seq;
}
