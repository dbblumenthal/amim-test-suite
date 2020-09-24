// Various probabilistic classes and functions, mostly about permutations
// should implement permutations as a class

#ifndef __RANDOM_H
#define __RANDOM_H

#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

inline int urand(int n) { // uniform random between 0 and n-1
  return (int) ((rand() / (((double) RAND_MAX) + 1)) * n);
}

void idPermutation(int *perm, int n); // set perm to identity
void randomizePermutation(int *perm, int n); // scramble perm
void printPermutation(ostream& str, int *perm, int n);

/*
 * PermutationGenerator class
 * yields a sequence of random permutations
 * first permutation returned is guaranteed to be Id
 * should implement version that enumerates all permutations
 */

class PermutationGenerator {
  int _count; // number of permutations already generated
  int _nPerms; // total number of permutations to generate
  bool _verbose; // if true, write . to cerr every _dotPeriod permutations
  int _dotPeriod;
  virtual void updatePermutation() {
    randomizePermutation(_perm, _n);
  }
 protected:
  int _n; // permutation size
  int *_perm; // holds the current permutation
 public:
  PermutationGenerator(int n, int nPerms, bool verbose = false, 
		       int dotPeriod = 100)
    : _n(n), _nPerms(nPerms), _count(0), _verbose(verbose), 
    _dotPeriod(dotPeriod),_perm(new int[_n]) {}
  ~PermutationGenerator() { delete _perm; }
  int nPerms(void) { return _nPerms; }
  int size(void) { return _n; }
  int *get() { // main method, returns next permutation
    if (++_count > _nPerms) { // done generating
      if (_verbose) cerr << endl; 
      return NULL;
    }
    else {
      if (_verbose && _count % _dotPeriod == 0) {
	if (_count % (10 * _dotPeriod) == 0)
	  cerr << _count / (10 * _dotPeriod);
	else 
	  cerr << "."; // print dot
      }
      if (_count == 1)
	idPermutation(_perm, _n);
      else
	updatePermutation();
      return _perm;
    }
  }
};

// derived class: yields balanced permutations

class BalancedPermutationGenerator : public PermutationGenerator {
  int *_permS; // starting permutation, based on phenotype
  int _n0; // number of zeroes in phenotype
  int _n1; // number of ones
  int _balance; // number of zeroes that become ones
public:
  BalancedPermutationGenerator(int n, int nPerms, int *pheno, int balance = 0,
			       bool verbose = false, int dotPeriod = 100)
    : PermutationGenerator(n, nPerms, verbose, dotPeriod),
    _permS(new int[n]), _n0(0), _n1(0), _balance(balance) {

    for (int i = 0; i < n; i++) { // prepare starting permutation
      if (pheno[i] == 0)
	_permS[_n0++] = i;
      else if (pheno[i] == 1)
	_permS[n - ++_n1] = i;
    }
    if (_balance <= 0) // set to something close to optimal
      _balance = (_n0 < _n1 ? _n0 / 2 : _n1 / 2);
    if (_balance > _n0) _balance = _n0;
    if (_balance > _n1) _balance = _n1;
  }
  ~BalancedPermutationGenerator() { delete _permS; }

  void updatePermutation() {
    idPermutation(_perm, _n);
    randomizePermutation(_permS, _n0);
    randomizePermutation(_permS + _n - _n1, _n1);
    for (int i = 0; i < _balance; i++) {
      _perm[_permS[i]] = _permS[_n - i - 1];
      _perm[_permS[_n - i - 1]] = _permS[i];
    }
  }
};

// derived class: yields invariant permutations
// current algorithm only works if types are sorted

class InvariantPermutationGenerator : public PermutationGenerator {
  const vector<int> types; // types are integers from 0 to nTypes - 1
  const int nTypes;
  vector<int> typeCount;

 public:
  InvariantPermutationGenerator(int n, int nPerms, const vector<int>& ty, 
			  int nTy, bool verbose = false, int dotPeriod = 100)
    : PermutationGenerator(n, nPerms, verbose, dotPeriod),
    types(ty), nTypes(nTy), typeCount(nTypes) {

    for (int i = 0; i < n; i++)
      typeCount[types[i]]++;
  }
  void updatePermutation() {
    int i, sum = 0;
    idPermutation(_perm, _n);
    for (i = 0; i < nTypes; i++) {
      randomizePermutation(_perm + sum, typeCount[i]);
      sum += typeCount[i];
    }
  }
};

#endif
