#include "Data.h"
#include "cdf.h"

int FastDataSet::pheno0 = 0;
int FastDataSet::pheno1 = 1;

double FastDataSet::sd(void) {
  if (n <= 1) return 0;
  double mu = sumx / n;
  double var = sumxx / n - mu * mu;
  return var > 0 ? sqrt(var) : 0;
}

// F statistic for a vector with associated phenotypes
// phenotypes are required to be integers btwn 0 and nPheno - 1
double fstat(double *x, const vector<int> pheno, const int nPheno)
{
  vector<FastDataSet> vv(nPheno); // stores stats for each phenotype
  vector<int>::const_iterator pp;
  for (pp = pheno.begin(); pp != pheno.end(); pp++, x++)
    vv[*pp].insert(*x);
  double sum = 0, sumsq = 0, sumrss = 0;
  int i;
  for (i = 0; i < nPheno; i++) {
    sum += vv[i].sumx;
    sumsq += vv[i].sumxx;
    if (vv[i].n > 1)
      sumrss += (vv[i].sumxx - vv[i].sumx * vv[i].sumx / vv[i].n);
  }
  double rssnull = sumsq - sum * sum / pheno.size();
  double rssdiff = rssnull - sumrss;
  if (sumrss > 0)
    return rssdiff / (nPheno-1) / sumrss * (pheno.size() - nPheno);
  else
    return 0; // not quite correct, but will prevent overflow
}
