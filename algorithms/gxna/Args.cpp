#include <sstream>
#include <fstream>
#include "Args.h"

Args::Args() :
  namesFile("geneID2name"),
  edgeFile("human.gra"),
  mapFile(""),
  expFile(""),
  phenoFile(""),
  typesFile(""),
  name(""),
  version("000"),
  outputDir("."),
  dotPeriod(100),
  graphCount(25),
  maxOverlap(0.75),
  maxRows(250),
  runDOT(false),
  algoType(0),
  scoreType(-1),
  compare("all"),
  takeAbs(false),
  maxTscaled(false),
  flexSize(false),
  minSD(0),
  minDegree(0),
  depth(15),
  iterations(0),
  kMet(1),
  radius(0),
  scalingExponent(0.6),
  seed(5),
  shrink(false),
  shrinkDF(0),
  shrinkVar(1),
  nPerms(100),
  useInvariantPerms(false),
  balance(0)
  {}

int Args::parse(int argc, char *argv[]) {
  int n = 1;
  while (n <= argc - 2) {
    istringstream s(argv[n+1]);
    if (string(argv[n]) == "-argFile") {
      ifstream fs(argv[n+1]);
      fs >> *this;
    }
    else if (string(argv[n]) == "-namesFile")
      s >> namesFile;
    else if (string(argv[n]) == "-edgeFile")
      s >> edgeFile;
    else if (string(argv[n]) == "-mapFile")
      s >> mapFile;
    else if (string(argv[n]) == "-expFile")
      s >> expFile;
    else if (string(argv[n]) == "-phenoFile")
      s >> phenoFile;
    else if (string(argv[n]) == "-typesFile")
      s >> typesFile;
    else if (string(argv[n]) == "-name") {
      s >> name;
      expFile = name + ".exp";
      phenoFile = name + ".phe";
      typesFile = name + ".typ";
    }
    else if (string(argv[n]) == "-version")
      s >> version;
    else if (string(argv[n]) == "-outputDir")
      s >> outputDir;
    else if (string(argv[n]) == "-dotPeriod")
      s >> dotPeriod;
    else if (string(argv[n]) == "-graphCount")
      s >> graphCount;
    else if (string(argv[n]) == "-maxOverlap")
      s >> maxOverlap;
    else if (string(argv[n]) == "-maxRows")
      s >> maxRows;
    else if (string(argv[n]) == "-runDOT")
      s >> runDOT;
    else if (string(argv[n]) == "-algoType") {
      s >> algoType;
      if (algoType > 0 && minDegree == 0) minDegree = 1;
    }
    else if (string(argv[n]) == "-scoreType")
      s >> scoreType;
    else if (string(argv[n]) == "-compare")
      s >> compare;
    else if (string(argv[n]) == "-takeAbs")
      s >> takeAbs;
    else if (string(argv[n]) == "-maxTscaled")
      s >> maxTscaled;
    else if (string(argv[n]) == "-flexSize")
      s >> flexSize;
    else if (string(argv[n]) == "-minSD")
      s >> minSD;
    else if (string(argv[n]) == "-minDegree")
      s >> minDegree;
    else if (string(argv[n]) == "-depth")
      s >> depth;
    else if (string(argv[n]) == "-iterations")
      s >> iterations;
    else if (string(argv[n]) == "-kMet")
      s >> kMet;
    else if (string(argv[n]) == "-radius") {
      s >> radius;
      if (radius > 0 && minDegree == 0) minDegree = 1;
    }
    else if (string(argv[n]) == "-scalingExponent")
      s >> scalingExponent;
    else if (string(argv[n]) == "-seed")
      s >> seed;
    else if (string(argv[n]) == "-shrink")
      s >> shrink;
    else if (string(argv[n]) == "-shrinkDF")
      s >> shrinkDF;
    else if (string(argv[n]) == "-shrinkVar")
      s >> shrinkVar;
    else if (string(argv[n]) == "-nPerms")
      s >> nPerms;
    else if (string(argv[n]) == "-useInvariantPerms")
      s >> useInvariantPerms;
    else if (string(argv[n]) == "-balance")
      s >> balance;
    else {
      cerr << "Args::parse: skipping unknown token " << argv[n] << endl;
    }
    n += 2;
  }
  return n;
}

istream& operator>>(istream& fs, Args& args) {
  string keyword;
  while (fs >> keyword) {
    if (keyword == "namesFile")
      fs >> args.namesFile;
    else if (keyword == "edgeFile")
      fs >> args.edgeFile;
    else if (keyword == "mapFile")
      fs >> args.mapFile;
    else if (keyword == "expFile")
      fs >> args.expFile;
    else if (keyword == "phenoFile")
      fs >> args.phenoFile;
    else if (keyword == "typesFile")
      fs >> args.typesFile;
    else if (keyword == "name") {
      fs >> args.name;
      args.expFile = args.name + ".exp";
      args.phenoFile = args.name + ".phe";
      args.typesFile = args.name + ".typ";
    }
    else if (keyword == "version")
      fs >> args.version;
    else if (keyword == "outputDir")
      fs >> args.outputDir;
    else if (keyword == "dotPeriod")
      fs >> args.dotPeriod;
    else if (keyword == "graphCount")
      fs >> args.graphCount;
    else if (keyword == "maxOverlap")
      fs >> args.maxOverlap;
    else if (keyword == "maxRows")
      fs >> args.maxRows;
    else if (keyword == "runDOT")
      fs >> args.runDOT;
    else if (keyword == "algoType") {
      fs >> args.algoType;
      if (args.algoType > 0 && args.minDegree == 0) args.minDegree = 1;
    }
    else if (keyword == "scoreType")
      fs >> args.scoreType;
    else if (keyword == "compare")
      fs >> args.compare;
    else if (keyword == "takeAbs")
      fs >> args.takeAbs;
    else if (keyword == "maxTscaled")
      fs >> args.maxTscaled;
    else if (keyword == "flexSize")
      fs >> args.flexSize;
    else if (keyword == "minSD")
      fs >> args.minSD;
    else if (keyword == "minDegree")
      fs >> args.minDegree;
    else if (keyword == "depth")
      fs >> args.depth;
    else if (keyword == "iterations")
      fs >> args.iterations;
    else if (keyword == "kMet")
      fs >> args.kMet;
    else if (keyword == "radius") {
      fs >> args.radius;
      if (args.radius > 0 && args.minDegree == 0) args.minDegree = 1;
    }
    else if (keyword == "scalingExponent")
      fs >> args.scalingExponent;
    else if (keyword == "seed")
      fs >> args.seed;
    else if (keyword == "shrink")
      fs >> args.shrink;
    else if (keyword == "shrinkDF")
      fs >> args.shrinkDF;
    else if (keyword == "shrinkVar")
      fs >> args.shrinkVar;
    else if (keyword == "nPerms")
      fs >> args.nPerms;
    else if (keyword == "useInvariantPerms")
      fs >> args.useInvariantPerms;
    else if (keyword == "balance")
      fs >> args.balance;
  }
  return fs;
}

ostream& operator<<(ostream& log, Args& args) {
  log << "namesFile " << args.namesFile << endl;
  log << "edgeFile " << args.edgeFile << endl;
  log << "mapFile " << args.mapFile << endl;
  log << "expFile " << args.expFile << endl;
  log << "phenoFile " << args.phenoFile << endl;
  log << "typesFile " << args.typesFile << endl;
  log << "name " << args.name << endl;
  log << "version " << args.version << endl;
  log << "outputDir " << args.outputDir << endl;
  log << "dotPeriod " << args.dotPeriod << endl;
  log << "graphCount " << args.graphCount << endl;
  log << "maxOverlap " << args.maxOverlap << endl;
  log << "maxRows " << args.maxRows << endl;
  log << "runDOT " << args.runDOT << endl;
  log << "algoType " << args.algoType << endl;
  log << "scoreType " << args.scoreType << endl;
  log << "compare " << args.compare << endl;
  log << "takeAbs " << args.takeAbs << endl;
  log << "maxTscaled " << args.maxTscaled << endl;
  log << "flexSize " << args.flexSize << endl;
  log << "minSD " << args.minSD << endl;
  log << "minDegree " << args.minDegree << endl;
  log << "depth " << args.depth << endl;
  log << "iterations " << args.iterations << endl;
  log << "kMet " << args.kMet << endl;
  log << "radius " << args.radius << endl;
  log << "scalingExponent " << args.scalingExponent << endl;
  log << "seed " << args.seed << endl;
  log << "shrink " << args.shrink << endl;
  log << "shrinkDF " << args.shrinkDF << endl;
  log << "shrinkVar " << args.shrinkVar << endl;
  log << "nPerms " << args.nPerms << endl;
  log << "useInvariantPerms " << args.useInvariantPerms << endl;
  log << "balance " << args.balance << endl;
  return log;
}
