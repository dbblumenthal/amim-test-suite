#ifndef __ARGS_H
#define __ARGS_H

#include <iostream>
#include <string>
using namespace std;

class Args {
 public:
  string namesFile;
  string edgeFile;
  string mapFile;
  string expFile;
  string phenoFile;
  string typesFile;
  string name;
  string version;
  string outputDir;
  int dotPeriod;
  int graphCount;
  double maxOverlap;
  int maxRows;
  bool runDOT;
  int algoType;
  int scoreType;
  string compare;
  bool takeAbs;
  bool maxTscaled;
  bool flexSize;
  double minSD;
  int minDegree;
  int depth;
  int iterations;
  double kMet;
  int radius;
  double scalingExponent;
  int seed;
  bool shrink;
  double shrinkDF;
  double shrinkVar;
  int nPerms;
  bool useInvariantPerms;
  int balance;

  Args();
  int parse(int argc, char *argv[]);
};

istream& operator>>(istream&, Args&);
ostream& operator<<(ostream&, Args&);

#endif
