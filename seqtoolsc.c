#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

//             ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz
#define TRANS "TVGHEFCDIJMLKNOPQYWAABSXRZ[\\]^_`tvghefcdijmlknopqywaabsxrz"
#define TRANS_OFFSET 65
#define TRANS_LEN 57

char* get_revcomp(char *input);
char get_char_comp(char c);
int *get_diffs_simple(char *cons, char *seqs[], int n_seqs);
double *get_diffs_frac_simple(char *cons, char *seqs[], int n_seqs);
double **get_diffs_frac_binned(char *cons, char *seqs[], int n_seqs, int seq_len, int bins);
char *transfer_gaps(char *gapped_seq, char *inseq, char gap_char1, char gap_char2);
char **transfer_gaps_multi(int n_seqs, char *gapped_seqs[], char *inseqs[], char gap_char1, char gap_char2);


// Return the reverse complement of a sequence.
// Makes a new copy of the string, so the original is not modified.
char* get_revcomp(char *input) {
  int length = strlen(input);
  char *output = malloc(sizeof(char) * length + 1);
  int i, j;
  for (i = 0, j = length - 1; i < length && j >= 0; i++, j--) {
    output[j] = get_char_comp(input[i]);
  }
  output[length] = '\0';
  return output;
}


// Return the complement of a base.
// Uses a simple lookup table: a string with the complements of all possible sequence characters.
char get_char_comp(char c) {
  int i = c - TRANS_OFFSET;
  if (i < 0 || i > TRANS_LEN) {
    return c;
  } else {
    return TRANS[i];
  }
}


/* Take an existing alignment and consensus and compute the number of differences between each
 * sequence and the consensus.
 * Known bugs:
 * 1. Counts no differences in the following sequences:
 *   consensus: GA---CA
 *   seq 1:     GA----A
 *   seq 2:     GA--ACA
 * 2. If a sequence starts with a gap, each base in the gap will be counted as a diff.
 */
int *get_diffs_simple(char *cons, char *seqs[], int n_seqs) {
  int *diffs = malloc(sizeof(int) * n_seqs);
  int i = 0;
  // Uppercase the consensus.
  while (cons[i] != 0) {
    cons[i] = toupper(cons[i]);
    i++;
  }
  // Loop through the sequences in the alignment.
  for (i = 0; i < n_seqs; i++) {
    int in_gap;
    diffs[i] = 0;
    int j = 0;
    // Compare each base of the sequence to the consensus.
    while (seqs[i][j] != 0 && cons[j] != 0) {
      if (cons[j] != '-' && seqs[i][j] != '-') {
        in_gap = 0;
      }
      if (toupper(seqs[i][j]) != cons[j]) {
        if (!in_gap) {
          diffs[i]++;
        }
      }
      if (cons[j] == '-' || seqs[i][j] == '-') {
        in_gap = 1;
      }
      j++;
    }
  }
  return diffs;
}


// Convert the output of get_diffs_simple() from raw diff counts to fractions of the total sequence
// lengths.
//TODO: Don't count gaps in sequence length.
double *get_diffs_frac_simple(char *cons, char *seqs[], int n_seqs) {
  int *diffs = get_diffs_simple(cons, seqs, n_seqs);
  double *fracs = malloc(sizeof(double) * n_seqs);
  int i;
  for (i = 0; i < n_seqs; i++) {
    int j = 0;
    while (seqs[i][j] != 0 && cons[j] != 0) {
      j++;
    }
    fracs[i] = (double)diffs[i]/j;
  }
  return fracs;
}


/* Take an existing alignment and consensus and compute the number of differences between each
 * sequence and the consensus. Break each sequence into bins and tally the differences in each bin.
 * Known bugs:
 * 1. counts no differences in the following sequences:
 *   consensus: GA---CA
 *   seq 1:     GA----A
 *   seq 2:     GA--ACA
 * 2. If a bin starts with a gap, each base in the gap will be counted as a diff.
 */
int **get_diffs_binned(char *cons, char *seqs[], int n_seqs, int seq_len, int bins) {
  int bin_size = (int)round((float)seq_len/bins);
  // Initialize the diffs 2d array.
  int **diffs = malloc(sizeof(int*) * n_seqs);
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    diffs[i] = malloc(bins * sizeof(int));
    for (j = 0; j < bins; j++) {
      diffs[i][j] = 0;
    }
  }
  // Uppercase the consensus.
  while (cons[i] != 0) {
    cons[i] = toupper(cons[i]);
    i++;
  }
  int bin, in_gap;
  // Loop through the sequences in the alignment.
  for (i = 0; i < n_seqs; i++) {
    j = 0;
    // Compare each base of the sequence to the consensus.
    while (seqs[i][j] != 0 && cons[j] != 0) {
      bin = j/bin_size;
      if (bin >= bins) {
        break;
      }
      if (cons[j] != '-' && seqs[i][j] != '-') {
        in_gap = 0;
      }
      if (toupper(seqs[i][j]) != cons[j]) {
        if (!in_gap) {
          diffs[i][bin]++;
        }
      }
      if (cons[j] == '-' || seqs[i][j] == '-') {
        in_gap = 1;
      }
      j++;
    }
  }
  return diffs;
}


// Convert the output of get_diffs_binned() from raw diff counts to fractions of the total bin
// lengths.
//TODO: Don't count gaps in bin length.
double **get_diffs_frac_binned(char *cons, char *seqs[], int n_seqs, int seq_len, int bins) {
  int bin_size = (int)round((float)seq_len/bins);
  int **diffs = get_diffs_binned(cons, seqs, n_seqs, seq_len, bins);
  double **fracs = malloc(sizeof(double*) * n_seqs);
  int i;
  for (i = 0; i < n_seqs; i++) {
    fracs[i] = malloc(sizeof(double) * bins);
    // Create and init array of lengths of the bins.
    int bin_lengths[bins];
    int bin;
    for (bin = 0; bin < bins; bin++) {
      bin_lengths[bin] = 0;
    }
    // Tally size of each bin.
    int j = 0;
    while (seqs[i][j] != 0 && cons[j] != 0) {
      int bin = j/bin_size;
      if (bin >= bins) {
        break;
      }
      bin_lengths[bin]++;
      j++;
    }
    // For each bin, calculate the diff fraction = diffs / bin_length.
    for (bin = 0; bin < bins; bin++) {
      fracs[i][bin] = (double)diffs[i][bin]/bin_lengths[bin];
      // printf("bin %d: %d / %d = %f\t", bin, diffs[i][bin], bin_lengths[bin], fracs[i][bin]);
    }
    // printf("\n");
  }
  return fracs;
}


// Take an input sequence and insert gaps according to another, already-aligned sequence with gaps.
// Input strings must be null-terminated. "gap_char1" is the character used for gaps in
// "gapped_seq", and "gap_char2" is the gap character in "inseq".
// N.B.: The ungapped length of "gapped_seq" must be equal to the length of "inseq".
char *transfer_gaps(char *gapped_seq, char *inseq, char gap_char1, char gap_char2) {
  if (gap_char1 == 0) {
    gap_char1 = '-';
  }
  if (gap_char2 == 0) {
    gap_char2 = '-';
  }
  int gapped_len = strlen(gapped_seq);
  char *outseq = malloc(sizeof(char) * gapped_len + 1);

  // Transfer characters from inseq to outseq, except when gapped_seq has a gap at that spot
  // (insert a gap there instead).
  int g, o, i;
  for (g = 0, o = 0, i = 0; g < gapped_len; g++, o++) {
    if (gapped_seq[g] == gap_char1) {
      outseq[o] = gap_char2;
    } else {
      outseq[o] = inseq[i];
      i++;
    }
  }
  outseq[gapped_len] = '\0';

  return outseq;
}


// Wrapper for transfer_gaps() when operating on a set of sequences at once.
char **transfer_gaps_multi(int n_seqs, char *gapped_seqs[], char *inseqs[], char gap_char1,
                           char gap_char2) {
  char **outseqs = malloc(sizeof(char *) * n_seqs);
  int i;
  for (i = 0; i < n_seqs; i++) {
    outseqs[i] = transfer_gaps(gapped_seqs[i], inseqs[i], gap_char1, gap_char2);
  }
  return outseqs;
}
