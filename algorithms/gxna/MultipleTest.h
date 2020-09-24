#ifndef __MULTIPLETEST_H
#define __MULTIPLETEST_H

/*
 * The class Stat must have:
 *   a method newPermutation() to be called for every new permutation
 *   an operator(): Stat(object) to compute the score/statistic for an object
 *     (and for the permutation last processed by newPermutation()
 *   a print(ostream&, Object) method for output
 *
 * For two-sided testing (typical), operator() should
 *   return the absolute value of the test statistic
 *
 * The first permutation generated MUST be id
 *
 * verboseLevel can be used for debugging:
 *  0 = completely quiet
 *  1 = report testing has started
 *  2 = in addition, print permutation and top object
 */

#include "Random.h"
#include "Data.h"
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
using namespace std;

class CompareObjects {
  const vector<double> _origT;
 public:
  CompareObjects(const vector<double>& origT) : _origT(origT) {}
  bool operator() (const int i, const int j) const {
    return _origT[i] > _origT[j];
  }
};

template<class Object, class Stat> class MultipleTest {
  Stat *stat;
  int nObjects;
  int verboseLevel;
  typename vector<Object>::const_iterator objects;
  vector<double> rawP;
  vector<double> adjP;
  vector<double> origT;
  vector<int> rank;

 public:
  int size() const { return nObjects; }
  double getRawP(int n) const { return rawP[n]; }
  double getAdjP(int n) const { return adjP[n]; }
  Object getObject(int n) const { return objects[rank[n]]; } // n-th ranked object
  MultipleTest(Stat& st, const vector<Object>& objs) : 
    stat(&st),
    nObjects(objs.size()),
    verboseLevel(1),
    objects(objs.begin()),
    rawP(nObjects),
    adjP(nObjects),
    origT(nObjects),
    rank(nObjects)
    {
      for (int i = 0; i < nObjects; i++) {
	rawP[i] = adjP[i] = 1; // since first permutation is id
	rank[i] = i;
      }
    }

  ostream& print(ostream& str) {
    typename multimap<double,Object>::iterator mi;
    for (int i = 0; i < nObjects; i++) { // go thru all objects
      str << i << ' ';
      stat->print(str, objects[rank[i]]);
      // str << ' ' << origT[rank[i]];
      // does not print test statistics because it does not know its sign
      // also it may not be the appropriate thing to print i.e.
      //   if testStat = T - meanT, we may want to print T instead of testStat
      // if desired, stat->print should print it
      str << ' ' << rawP[i] << ' ' << adjP[i] << endl;
    }
    return str;
  }

  void maxT(PermutationGenerator& pg) {
    if (verboseLevel >= 1)
      cerr << "maxT: testing " << nObjects << " objects\n";
    int nPerms = 0, i, argMax;
    double currentMax;
    while (int *currentPerm = pg.get()) { // generate new permutation
      ++nPerms;
      stat->newPermutation(currentPerm);
      if (nPerms == 1) { // first perm housekeeping
	for (i = 0; i < nObjects; i++)
	  origT[i] = (*stat)(objects[i]);
	sort(rank.begin(), rank.end(), CompareObjects(origT));
	argMax = 0;
	currentMax = origT[rank[argMax]];
      }
      else {
	currentMax = -(numeric_limits<double>::max());
	argMax = -1;
	for (i = nObjects - 1; i >= 0; i--) { // go thru all objects
	  double T = (*stat)(objects[rank[i]]);
	  if (T >= origT[rank[i]]) // update raw count
	    rawP[i]++;
	  if (currentMax < T) {
	    currentMax = T;
	    argMax = i;
	  }
	  if (currentMax >= origT[rank[i]]) // update adjusted count
	    adjP[i]++;
	}
      }
      if (verboseLevel >= 2) {
	cerr << "maxT: " << nPerms << ' ' << argMax << ' ';
	stat->print(cerr, objects[rank[argMax]]);
	cerr << ' ' << currentMax << ' ';
	printPermutation(cerr, currentPerm, pg.size());
	cerr << endl;
      }
    } // done generating permutations
    for (i = 0; i < nObjects; i++) { // compute p-values
      rawP[i] /= nPerms;
      adjP[i] /= nPerms;
      if (i > 0 && adjP[i] < adjP[i-1]) // enforce monotonicity
        adjP[i] = adjP[i-1];
    }
  } // end of maxT

  void maxTscaled(PermutationGenerator& pg) {
    bool PrintMuAndSigma = false;
    if (verboseLevel >= 1)
      cerr << "maxTscaled: testing " << nObjects << " objects\n";
    int nPerms = 0, i, j, tt = 0;
    vector<double> Tarray(nObjects * pg.nPerms());
    vector<FastDataSet> Tdata(nObjects);
    while (int *currentPerm = pg.get()) { // compute all T's first
      ++nPerms;
      stat->newPermutation(currentPerm);
      for (i = 0; i < nObjects; i++) {
	double T = (*stat)(objects[i]);
	Tarray[tt] = T;
	tt++;
	Tdata[i].insert(T);
      }
    }
    for (i = 0; i < nObjects; i++) { // rescale T's
      double mu = Tdata[i].mean();
      double sigma = Tdata[i].sd();
      if (PrintMuAndSigma) {
	stat->print(cout, objects[i]);
	cout << " " << mu << " " << sigma << endl;
      }
      origT[i] = Tarray[i] - mu;
      if (sigma > 0)
	origT[i] /= sigma;
      rank[i] = i;
      for (j = 0; j < nPerms; j++) {
	Tarray[i + j * nObjects] -= mu;
	if (sigma > 0)
	  Tarray[i + j * nObjects] /= sigma;
      }
    }
    sort(rank.begin(), rank.end(), CompareObjects(origT)); // rerank
    for (j = 1; j < nPerms; j++) { // compute p-values
      double currentMax = -(numeric_limits<double>::max());
      int argMax = -1;
      for (i = nObjects - 1; i >= 0; i--) { // go thru all objects
	double T = Tarray[rank[i] + j * nObjects];
	if (T >= origT[rank[i]]) // update raw count
	  rawP[i]++;
	if (currentMax < T) {
	  currentMax = T;
	  argMax = i;
	}
	if (currentMax >= origT[rank[i]]) // update adjusted count
	  adjP[i]++;
      }
      if (verboseLevel >= 2) {
	cerr << "maxTscaled: " << j << ' ' << argMax << ' ';
	stat->print(cerr, objects[rank[argMax]]);
	cerr << ' ' << currentMax << endl;
      }
    }
    for (i = 0; i < nObjects; i++) {
      rawP[i] /= nPerms;
      adjP[i] /= nPerms;
      if (i > 0 && adjP[i] < adjP[i-1]) // enforce monotonicity
        adjP[i] = adjP[i-1];
    }
  } // end of maxTscaled
};

template<class Object, class Stat> ostream& operator<<(ostream& str, 
	MultipleTest<Object,Stat>& mt) {
  return mt.print(str);
}

#endif
