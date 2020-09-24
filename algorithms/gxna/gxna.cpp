// main executable
// balance is currently not used

#include "GEData.h"
#include "Graph.h"
#include "Args.h"
#include "MultipleTest.h"

int main(int argc, char *argv[]) {
  cerr << "GXNA Version 2.0 (Beta) released March 25, 2008\n";
  Args arg;
  ifstream argfs("defaults.txt");
  argfs >> arg;
  arg.parse(argc, argv);
  if (arg.mapFile == "") {
    cerr << "You must specify a mapFile that maps probes to Entrez id.\n";
    exit(1);
  }

  GEData ged(arg); // create main object
  srand(arg.seed); // seed random number generator
  string argsFileName = arg.outputDir + "/" + arg.name + "_"
    + arg.version + ".arg";
  ofstream fs0(argsFileName.c_str());
  fs0 << arg;

  /*
  PermutationGenerator *ppg;
  if (arg.balance == 0)
    ppg = &PermutationGenerator(ged.getNArrays(), arg.nPerms, true);
  else
    ppg = &BalancedPermutationGenerator(ged.getNArrays(), arg.nPerms, 
					ged.getPhenotype(), arg.balance, true);
  */

  // do multiple testing
  if (arg.useInvariantPerms) {
    InvariantPermutationGenerator ppg(ged.getNArrays(), arg.nPerms, 
				      ged.arrayType, ged.nTypes, true, arg.dotPeriod);
    if (arg.algoType == 0)
      ged.testBall(ppg);
    else
      ged.testSubgraphs(ppg);
  }
  else {
    PermutationGenerator ppg(ged.getNArrays(), arg.nPerms, true, arg.dotPeriod);
    if (arg.algoType == 0)
      ged.testBall(ppg);
    else
      ged.testSubgraphs(ppg);
  }
  return 0;
}
