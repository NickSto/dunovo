#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

// N.B. This defines the valid bases, but it's also effectively defined in the switches in
// get_votes_simple(), get_votes_qual(), and get_base_prime(), and in the constant IUPAC_BASES.
#define N_BASES 6
const char *BASES = "ACGTN-";
/* A  C   G   T   N   -     A: 2    Compute IUPAC ambiguous base character by representing each base
A  4  6  10  14  22  26     C: 3    with a prime and multiplying. Then use a lookup table (an array
C     9  15  21  33  39     G: 5    where the index is the product of the two primes).
G        25  35  55  65     T: 7
T            49  77  91     N: 11
N               121 143     -: 13    1         2         3         4         5         6         7
-                   169    01234567890123456789012345678901234567890123456789012345678901234567890*/
const char *IUPAC_BASES = "N...A.M..CR...WS.....YN..GN......N.K...N.........T.....N.........N....."
//                                   8         9        10        11        12        13        14
                           "......N.............N.............................N..................."
//                                  15        16        17
                           "..N.........................-";
#define THRES_DEFAULT 0.5
#define WIN_LEN 4

int **get_votes_simple(char *align[], int n_seqs, int seq_len);
int **get_votes_qual(char *align[], char *quals[], int n_seqs, int seq_len, char thres);
int **init_votes(int seq_len);
void free_votes(int *votes[], int seq_len);
void print_votes(char *consensus, int *votes[], int seq_len);
char *rm_gaps(char *consensus, int cons_len);
char *build_consensus(int *votes[], int seq_len, double thres);
char *build_consensus_duplex(int *votes1[], int *votes2[], int seq_len, double thres);
char *build_consensus_duplex_simple(char *cons1, char *cons2, int gapped);
int get_base_prime(char base);
char *get_consensus(char *align[], char *quals[], int n_seqs, int seq_len, double thres,
                    char qual_thres, int gapped);
char *get_consensus_duplex(char *align1[], char *align2[], char *quals1[], char *quals2[],
                           int n_seqs1, int n_seqs2, int seq_len, double cons_thres,
                           char qual_thres, int gapped, char *method);


// Tally the different bases at each position in an alignment.
// Returns an array of arrays: for each position in the alignment, an array of the number of times
// each base occurs at that position. The order of bases is as in the "BASES" constant.
int **get_votes_simple(char *align[], int n_seqs, int seq_len) {
  int **votes = init_votes(seq_len);

  // Tally votes for each base.
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    for (j = 0; j < seq_len; j++) {
      // N.B.: Could write this without hardcoded literals, but it's about 40% slower.
      switch (toupper(align[i][j])) {
        case 'A':
          votes[j][0]++;
          break;
        case 'C':
          votes[j][1]++;
          break;
        case 'G':
          votes[j][2]++;
          break;
        case 'T':
          votes[j][3]++;
          break;
        case 'N':
          votes[j][4]++;
          break;
        case '-':
          votes[j][5]++;
          break;
      }
    }
  }

  return votes;
}


int **get_votes_qual(char *align[], char *quals[], int n_seqs, int seq_len, char thres) {
  int **votes = init_votes(seq_len);

  // Tally votes for each base.
  int i, j;
  for (i = 0; i < n_seqs; i++) {
    for (j = 0; j < seq_len; j++) {
      // Don't count bases whose quality is less than the threshold. Gaps are always counted.
      //TODO: What to use as the quality for gaps? Currently gaps are effectively always considered
      //      the highest quality, biasing toward gaps.
      if (align[i][j] != '-' && quals[i][j] < thres) {
        continue;
      }
      // N.B.: Could write this without hardcoded literals, but it's about 40% slower.
      switch (toupper(align[i][j])) {
        case 'A':
          votes[j][0]++;
          break;
        case 'C':
          votes[j][1]++;
          break;
        case 'G':
          votes[j][2]++;
          break;
        case 'T':
          votes[j][3]++;
          break;
        case 'N':
          votes[j][4]++;
          break;
        case '-':
          votes[j][5]++;
          break;
      }
    }
  }

  return votes;
}


int *init_gap_qual_window(char *quals, int seq_len, int *win_edge) {
  // This does the initial fill of the "window" array, adding the first WIN_LEN bases to the right
  // side.
  int *window = malloc(sizeof(int) * WIN_LEN * 2);
  // Fill left side with -1's (no quality information).
  int i;
  for (i = 0; i < WIN_LEN; i++) {
    window[i] = -1;
  }
  // Fill right side with first WIN_LEN quality scores (or -1 if seq_len < WIN_LEN).
  for (i = WIN_LEN; i < WIN_LEN*2; i++) {
    if (i-WIN_LEN < seq_len) {
      window[i] = quals[i-WIN_LEN];
    } else {
      window[i] = -1;
    }
  }
  (*win_edge) = WIN_LEN - 1;
  return window;
}


int get_gap_qual(int *window, int *win_edge, int pos, char *quals, int seq_len, char gap_char) {
  /* Implementation:
   * "window" is an array of length 2*WIN_LEN, holding the quality scores of the bases WIN_LEN from
   * from the current one in both directions. At an edge, when there is no base to fill a given slot
   * in "window", -1 is a sentinel for "empty". For example, if we're at the 2nd base in the
   * sequence ("pos" 1), and WIN_LEN is 4, "window" should look something like this:
   *   base coordinates            0  2  3  4  5
   *   array values     [-1|-1|-1|24|32|29|36|33]
   * Each time the calling loop moves one base forward, call this function to add another quality
   * value to the right side of the array, and shift the rest of the values left. At the same time,
   * it will calculate a weighted average of the values, with the bases closest to the center (the
   * base of interest) weighed highest.
   */
  // Find the next quality score that's not a gap.
  char next_qual = gap_char;
  while (next_qual == gap_char) {
    (*win_edge)++;
    if ((*win_edge) < seq_len) {
      next_qual = quals[(*win_edge)];
    } else {
      (*win_edge) = seq_len;
      next_qual = -1;
    }
  }
  // Shift all the quality scores left and average them.
  //TODO: Do a weighted average.
  int total = 0;
  int scores = 0;
  int last_qual;
  int i;
  for (i = WIN_LEN*2 - 1; i >= 0; i--) {
    last_qual = window[i];
    window[i] = next_qual;
    next_qual = last_qual;
    if (window[i] != -1) {
      total += window[i];
      scores++;
    }
  }
  return total/scores;
}


void print_window(int *window) {
  printf("[");
  int i;
  for (i = 0; i < WIN_LEN*2; i++) {
    if (i == WIN_LEN*2 - 1) {
      printf("%2d]\n", window[i]);
    } else {
      printf("%2d|", window[i]);
    }
  }
}


int **init_votes(int seq_len) {
  int **votes = malloc(sizeof(int *) * seq_len);
  int i, j;
  for (i = 0; i < seq_len; i++) {
    votes[i] = malloc(sizeof(int) * N_BASES);
    for (j = 0; j < N_BASES; j++) {
      votes[i][j] = 0;
    }
  }
  return votes;
}


void free_votes(int *votes[], int seq_len) {
  int i;
  for (i = 0; i < seq_len; i++) {
    free(votes[i]);
  }
  free(votes);
}


void print_votes(char *consensus, int *votes[], int seq_len) {
  int i, j;
  printf("   ");
  for (j = 0; j < N_BASES; j++) {
    printf(" %c ", BASES[j]);
  }
  printf("\n");
  for (i = 0; i < seq_len; i++) {
    printf("%c: ", consensus[i]);
    for (j = 0; j < N_BASES; j++) {
      if (votes[i][j]) {
        printf("%2d ", votes[i][j]);
      } else {
        printf("   ");
      }
    }
    printf("\n");
  }
}


// Take a consensus sequence which may have gaps ('-' characters) and remove them to produce the
// actual final sequence. "cons_len" should be the length of the original, gapped, sequence.
char *rm_gaps(char *consensus, int cons_len) {
  char *output = malloc(sizeof(char) * cons_len + 1);
  int i;
  int j = 0;
  for (i = 0; i < cons_len; i++) {
    if (consensus[i] != '-') {
      output[j] = consensus[i];
      j++;
    }
  }
  output[j] = '\0';
  return output;
}


char *build_consensus(int *votes[], int seq_len, double thres) {
  char *consensus = malloc(sizeof(char) * seq_len + 1);
  
  int i, j;
  for (i = 0; i < seq_len; i++) {
    int total = 0;
    int max_vote = 0;
    char max_base = 'N';
    for (j = 0; j < N_BASES; j++) {
      total += votes[i][j];
      if (votes[i][j] > max_vote) {
        max_vote = votes[i][j];
        max_base = BASES[j];
      }
      if (total == 0) {
        consensus[i] = 'N';
      } else if ((double)max_vote/total > thres) {
        consensus[i] = max_base;
      } else {
        consensus[i] = 'N';
      }
    }
  }

  consensus[seq_len] = '\0';
  return consensus;
}


// Build a consensus sequence from two alignments by weighting each equally and considering only
// the frequency of each base in each alignment.
char *build_consensus_duplex(int *votes1[], int *votes2[], int seq_len, double thres) {
  char *consensus = malloc(sizeof(char) * seq_len + 1);

  int i, j;
  for (i = 0; i < seq_len; i++) {
    // Sum the total votes at this position.
    /*TODO: This does an extra loop through the votes to get the total so it can calculate actual
     *      frequencies in the second pass. Technically, this information could be gathered when
     *      originally tallying the votes in the get_votes functions. Or, the total could be assumed
     *      to be n_seqs if every base always contributes a vote (even when it's not in "ACGTN-").
     */
    int total1 = 0;
    for (j = 0; j < N_BASES; j++) {
      total1 += votes1[i][j];
    }
    int total2 = 0;
    for (j = 0; j < N_BASES; j++) {
      total2 += votes2[i][j];
    }
    double max_freq = 0.0;
    char max_base = 'N';
    for (j = 0; j < N_BASES; j++) {
      // Get the frequency of each base.
      double freq1;
      if (total1 > 0) {
        freq1 = (double)votes1[i][j]/total1;
      }
      double freq2;
      if (total2 > 0) {
        freq2 = (double)votes2[i][j]/total2;
      }
      // frequency of the base = average of frequencies in the two sequences
      double avg_freq;
      if (total1 == 0 && total2 == 0) {
        avg_freq = -1.0;
      } else if (total1 == 0) {
        avg_freq = freq2;
      } else if (total2 == 0) {
        avg_freq = freq1;
      } else {
        avg_freq = (freq1 + freq2) / 2;
      }
      // Track the highest frequency seen.
      if (avg_freq > max_freq) {
        max_freq = avg_freq;
        max_base = BASES[j];
      }
    }
    if (max_freq > thres) {
      consensus[i] = max_base;
    } else {
      consensus[i] = 'N';
    }
  }

  consensus[seq_len] = '\0';
  return consensus;
}


// "cons1" and "cons2" must be null-terminated strings of equal lengths.
char *build_consensus_duplex_simple(char *cons1, char *cons2, int gapped) {
  int seq_len = strlen(cons1);
  char *cons = malloc(sizeof(char) * seq_len + 1);
  int i = 0;
  int base_prime1, base_prime2;
  while (cons1[i] != '\0' && cons2[i] != '\0') {
    base_prime1 = get_base_prime(cons1[i]);
    base_prime2 = get_base_prime(cons2[i]);
    cons[i] = IUPAC_BASES[base_prime1*base_prime2];
    i++;
  }
  cons[seq_len] = '\0';
  if (gapped) {
    return cons;
  } else {
    return rm_gaps(cons, seq_len);
  }
}


int get_base_prime(char base) {
  switch (base) {
    case 'A':
      return 2;
    case 'C':
      return 3;
    case 'G':
      return 5;
    case 'T':
      return 7;
    case 'N':
      return 11;
    case '-':
      return 13;
    default:
      return 0;
  }
}


// Convenience function to create a consensus in one step.
// Give 0 as "quals" to not use quality scores, and -1.0 as "cons_thres" to use the default
// consensus threshold when evaluating base votes.
char *get_consensus(char *align[], char *quals[], int n_seqs, int seq_len, double cons_thres,
                    char qual_thres, int gapped) {
  if (cons_thres == -1.0) {
    cons_thres = THRES_DEFAULT;
  }
  int **votes;
  if (quals == 0) {
    votes = get_votes_simple(align, n_seqs, seq_len);
  } else {
    votes = get_votes_qual(align, quals, n_seqs, seq_len, qual_thres);
  }
  char *consensus_gapped = build_consensus(votes, seq_len, cons_thres);
  char *consensus;
  if (gapped) {
    consensus = consensus_gapped;
  } else {
    consensus = rm_gaps(consensus_gapped, seq_len);
  }
  free_votes(votes, seq_len);
  return consensus;
}


char *get_consensus_duplex(char *align1[], char *align2[], char *quals1[], char *quals2[],
                           int n_seqs1, int n_seqs2, int seq_len, double cons_thres,
                           char qual_thres, int gapped, char *method) {
  if (cons_thres == -1.0) {
    cons_thres = THRES_DEFAULT;
  }
  int **votes1;
  int **votes2;
  if (quals1 == 0 || quals2 == 0) {
    votes1 = get_votes_simple(align1, n_seqs1, seq_len);
    votes2 = get_votes_simple(align2, n_seqs2, seq_len);
  } else {
    votes1 = get_votes_qual(align1, quals1, n_seqs1, seq_len, qual_thres);
    votes2 = get_votes_qual(align2, quals2, n_seqs2, seq_len, qual_thres);
  }
  char *consensus_gapped;
  if (!strncmp(method, "freq", 4)) {
    consensus_gapped = build_consensus_duplex(votes1, votes2, seq_len, cons_thres);
  } else if (!strncmp(method, "iupac", 5)) {
    char *cons1 = build_consensus(votes1, seq_len, cons_thres);
    char *cons2 = build_consensus(votes2, seq_len, cons_thres);
    consensus_gapped = build_consensus_duplex_simple(cons1, cons2, 1);
  } else {
    return "";
  }
  char *consensus;
  if (gapped) {
    consensus = consensus_gapped;
  } else {
    consensus = rm_gaps(consensus_gapped, seq_len);
  }
  free_votes(votes1, seq_len);
  free_votes(votes2, seq_len);
  return consensus;
}


int main(int argc, char *argv[]) {
  char **align = malloc(sizeof(char *) * (argc-1));

  int seq_len = INT_MAX;
  int i;
  for (i = 1; i < argc; i++) {
    if (strlen(argv[i]) < seq_len) {
      seq_len = strlen(argv[i]);
    }
    align[i-1] = argv[i];
  }

  if (argc <= 1) {
    return 1;
  }

  char *quals = align[0];
  seq_len = strlen(quals);
  int *win_edge = malloc(sizeof(int));
  int *window = init_gap_qual_window(quals, seq_len, win_edge);
  print_window(window);
  int gap_qual = get_gap_qual(window, win_edge, 0, quals, seq_len, ' ');
  print_window(window);
  printf("%d\n", gap_qual);
  gap_qual = get_gap_qual(window, win_edge, 0, quals, seq_len, ' ');
  print_window(window);
  printf("%d\n", gap_qual);
  return 0;

  int **votes = get_votes_simple(align, argc-1, seq_len);
  char *consensus = build_consensus(votes, seq_len, THRES_DEFAULT);
  print_votes(consensus, votes, seq_len);
  printf("%s\n", consensus);
  free_votes(votes, seq_len);

  return 0;
}
