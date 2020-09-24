// class for gene expression data

// should replace M by vector
// should make arrayType, nTypes private (must change gxna.cpp)
// may want to make gbest, groupScore local vars in their functions

#ifndef __GEDATA_H
#define __GEDATA_H

#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include "Graph.h"
#include "Random.h"
#include "Data.h"
#include "Args.h"
#include "MultipleTest.h"
using namespace std;

enum { ALGO_BASIC = 0, ALGO_GXNA = 1 };
enum { SCORE_UNKNOWN = -1, SCORE_SIGMA_T = 0, SCORE_T_SIGMA = 1, 
       SCORE_SIGMA_F = 2, SCORE_F_SIGMA = 3 };

inline bool isT(int i) { return i <= 1; } // score is based on T statistic
inline bool isF(int i) { return i >= 2; } // score is based on T statistic
inline bool isAll(string s) { return s == "all" || s == "ALL"; }

class GEData {
  // Basics
  Args& arg; // program arguments
  int nGenes;
  int nArrays;
  Graph gg; // gene interaction graph
  double *M; // matrix of M values
  vector<double> foldChange; // used only when comparing two phenotypes
  vector<double> sd; // gene standard deviations
  vector<double> sfac; // shrinkage factors (set to 1 if no shrinkage)

  vector<int> origPhenotype; // real phenotype
  vector<int> phenotype; // permuted phenotype
  map<string,int> phenoMap; // maps phenotype names to integers
  vector<string> n2pheno; // maps integers to phenotype names
  int nPheno; // number of different phenotypes

  // Various maps
  // Genes are integers, [LocusLink] IDs are strings
  vector<string> geneName;
  vector<string> geneID; // LocusLink
  map<string,string> probe2id;
  vector<int> probesPerGene;
  map<string,int> id2gene; // used in reading the graph gg

  // Used for graph search
  vector< vector<int> > subGraphs;
  vector<double> subGraphScore; // used only for output
  vector<int> gbest;

  double *groupScore; // used to compute TSigma

  // Following vars are used for multiple testing
  // T, signT and meanT are updated by newPermutation(), so they should
  //   only be used by the multiple testing folks i.e. operators() and
  //   their underlings
  // Do not use them for output etc.
  int permCount;
  vector<double> T; // statistics (T or F) for the current perm
  vector<int> signT; // sign of T stat for the current perm
  double meanT; // mean of all T or |T| for the current perm

  // Used for DOT printing
  int DOTPrinted; // counts no of DOT files already written
  string baseFileName; // basename for output files
  string outputFileName;
  string htmlFileName;
  string htmlFrameFileName;
  ofstream frameStream;

  // Input functions
  int readMap(const char *file); // returns no of probes read
  int readPhenotype(const char *file); // returns no of phenotypes read
  int readExpression(const char *file); // returns no of good lines read

  // helper scoring functions
  double vectorScore(double *values, const vector<int>& pheno); 
  double origScore(int gene) {
    return vectorScore(M + gene * nArrays, origPhenotype) * sfac[gene];
  }
  double groupSum(const vector<int>&, int);

 public:
  vector<int> arrayType; // used to generate invariant permutations
  int nTypes;

  GEData(Args&);
  ~GEData() { 
    delete M, groupScore;
  }
  int getNArrays() { return nArrays; }
  string name(int i) { return geneName[i]; }
  string geneInfo(int gene) { return geneName[gene] + " " + geneID[gene]; }
  ostream& print(ostream& str) { // prints the GEData object
    for (int i = 0; i < nGenes; i++)
      str << i << ' ' << geneInfo(i) << ' ' << sd[i] 
	  << ' ' << origScore(i) << endl;
    return str;
  }

  // functions used for multiple testing
  // Subgraph is the actual ball, in the pre-defined case
  // int is the root gene, in the adapted case
  void newPermutation(int *perm);
  double scoreSubGraph(const SubGraph& genes); // "pre-defined" score
  double operator() (int gene); // "adapted" score

  // output functions
  ostream& print(ostream& str, int gene); // prints the corresponding subgraph
  void makeDOT(int gene, int nDOT);
  void printResults(const MultipleTest<int,GEData>&);
  string graphFileName(int, const string&);
  void beginHTML();
  void addRow(int, const string&, const string&, int, double, double, double);
  void endHTML();

  // master computation functions
  void testBall(PermutationGenerator&);
  void testSubgraphs(PermutationGenerator&);
};

#endif
