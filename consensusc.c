#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

// N.B. This defines the valid bases, but it's also effectively defined in the switch in
// get_votes_simple()
#define N_BASES 6
const char *BASES = "ACGTN-";
#define THRES_DEFAULT 0.5

int **get_votes_simple(char *align[], int n_seqs, int seq_len);
double **get_freqs(int *votes[], int seq_len);
int **init_votes(int seq_len);
void free_votes(int *votes[], int seq_len);
double **init_freqs(int seq_len);
void free_freqs(double *freqs[], int seq_len);
void print_votes(char *consensus, int *votes[], int seq_len);
char *rm_gaps(char *consensus, int cons_len);
char *build_consensus(int *votes[], int seq_len, double thres);
char *get_consensus(char *align[], char *quals[], int n_seqs, int seq_len, double thres,
                    char qual_thres, int gapped);


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


// Convert the counts in the "votes" array returned by get_votes_simple() into base frequencies by
// dividing by the total number of bases at each position.
double **get_freqs(int *votes[], int seq_len) {
  double **freqs = init_freqs(seq_len);

  // This passes through the data twice, once to get the totals and once to calculate the freqs.
  /*TODO: One of these passes could be eliminated by either getting the totals in get_votes_simple()
   *      (but that would require breaking the separation btwn these functions) or by always having
   *      the total be equal to the number of sequences (only not true currently when a base isn't
   *      one of "ACGTN-").
   */
  int i, j;
  for (i = 0; i < seq_len; i++) {
    int total = 0;
    for (j = 0; j < N_BASES; j++) {
      total += votes[i][j];
    }
    for (j = 0; j < N_BASES; j++) {
      freqs[i][j] = votes[i][j]/total;
    }
  }

  return freqs;
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


double **init_freqs(int seq_len) {
  double **freqs = malloc(sizeof(double *) * seq_len);
  int i;
  for (i = 0; i < seq_len; i++) {
    freqs[i] = malloc(sizeof(double) * N_BASES);
  }
  return freqs;
}


void free_freqs(double *freqs[], int seq_len) {
  int i;
  for (i = 0; i < seq_len; i++) {
    free(freqs[i]);
  }
  free(freqs);
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

  int **votes = get_votes_simple(align, argc-1, seq_len);
  char *consensus = build_consensus(votes, seq_len, THRES_DEFAULT);
  print_votes(consensus, votes, seq_len);
  printf("%s\n", consensus);
  free_votes(votes, seq_len);

  return 0;
}
