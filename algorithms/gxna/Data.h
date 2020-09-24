// class for various statistical functions
// inline functions declared in header because of compiler quirk
//   (should fix some day...)

#ifndef __DATA_H
#define __DATA_H

#include <vector>
#include <math.h>
using namespace std;

#define MAX(a,b) (((a)>(b))?(a):(b))

// FastDataSet: very fast, very basic
// starts empty, grows by insertion; only tracks mu and sigma
// note the var estimate is the MLE (biased)

class FastDataSet {
  int n;
  double sumx;
  double sumxx; // sum of squares

 public:
  static int pheno0; // for tstatPheno
  static int pheno1;
  FastDataSet() { n = 0; sumx = sumxx = 0; }
  void insert(double x) { sumx += x; sumxx += x * x; ++n; } 
  double mean() { return n > 0 ? sumx / n : 0; }
  double var() { return n > 1 ? sumxx / n - (sumx / n) * (sumx / n) : 0; }
  double sd();
  double se() { return n > 1 ? sd() / sqrt((double) n) : 0; }
  int size() { return n; }
  static void setPheno(int n0, int n1) { pheno0 = n0; pheno1 = n1; }
  friend double tstat(FastDataSet&, FastDataSet&); // unequal variances
  friend double tstateq(FastDataSet&, FastDataSet&); // equal variances
  friend double fstat(double *, const vector<int>, const int);
};

// T statistic between two data sets, different variances assumed
inline double tstat(FastDataSet& ds0, FastDataSet& ds1) {
  if (ds0.n < 2 || ds1.n < 2) // sets too small
    return 0;
  else {
    double mu0 = ds0.sumx / ds0.n;
    double mu1 = ds1.sumx / ds1.n;
    double vardiff = (ds1.sumxx / ds1.n - mu1 * mu1) / ds1.n
      + (ds0.sumxx / ds0.n - mu0 * mu0) / ds0.n;
    return vardiff > 0 ? (mu1 - mu0) / sqrt(vardiff) : 0;
  }
}

// T statistic between two data sets, equal variances assumed
inline double tstateq(FastDataSet& ds0, FastDataSet& ds1) {
  if (ds0.n < 1 || ds1.n < 1) // sets too small
    return 0;
  else {
    double mu0 = ds0.sumx / ds0.n;
    double mu1 = ds1.sumx / ds1.n;
    int n = ds0.n + ds1.n;
    double mu = (ds0.sumx + ds1.sumx) / n;
    double var = ((ds0.sumxx + ds1.sumxx) / n - mu * mu) / n;
    return var > 0 ? (mu1 - mu0) / sqrt(var) : 0;
  }
}


// T statistic for a vector with associated phenotypes
// phenotypes 0 and 1 define the two data sets, everything else is ignored
inline double tstatPheno(double *x, const vector<int>& pheno) {
  FastDataSet ds0, ds1;
  vector<int>::const_iterator pp;
  for (pp = pheno.begin(); pp != pheno.end(); pp++, x++) {
    if (*pp == FastDataSet::pheno0)
      ds0.insert(*x);
    else if (*pp == FastDataSet::pheno1)
      ds1.insert(*x);
  }
  return tstat(ds0, ds1);
}

double fstat(double *, const vector<int>, const int);

#endif
