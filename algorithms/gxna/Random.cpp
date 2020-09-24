#include "Random.h"

// create ID permutation
void idPermutation(int *perm, int n) {
  for (int i = 0; i < n; i++)
    perm[i] = i;
}

// uniform random scramble
void randomizePermutation(int *perm, int n) {
  for (int i = 1; i < n; i++) {
    int j = urand(i+1); // uniform among 0, 1, ..., i
    int temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
  }
}

void printPermutation(ostream& str, int *perm, int n) {
  for (int i = 0; i < n; i++)
    str << perm[i] << ' ';
  str << endl;
}
