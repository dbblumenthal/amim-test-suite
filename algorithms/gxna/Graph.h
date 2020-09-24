#ifndef __GRAPH_H
#define __GRAPH_H

// Graph class
// Designed to coexist with GEData and Args, not to be self-standing
// should fix this at some point...
// (should probably make findSubgraph friend of GEData and Graph)

#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <sstream>
#include "Random.h"
#include "Args.h"
using namespace std;

typedef vector<int> SubGraph; // list of nodes; kind of a hack

// nodes are integers (0, 1, ...) and may have labels
// edges are represented indirectly:
//   each node has a set of its neighbors
// there is no constructor; must read graph from a file before using it
//   (this is a hack because GEData creates gg before it knows
//    the gene2ID map; should rewrite this)
// the caller of read() decides how many nodes there will be
// then edges are read by read() from the edgeFile
//
// currently undirected; will need some rewriting to make it directed
//

class Graph {
  int nNodes;
  int nEdges;
  vector<int> _degree; // node degrees
  vector<string> _label; // node labels
  vector< vector<int> > neighbors; // for each vertex a set of neighbors
  map< pair<int,int>, string > edgeType;
  map< pair<int,int>, string > edgeSource;
  void addEdge(int, int, string, string);
  
  // the following variables are set by setScores and used by findSubgraph
  vector<double> score; // a score for each node
  double meanScore;
  double sdScore;
  // compute score of a subgraph, given the sum of the scores of its vertices
  double scalingExponent;
  double sum2score(double sum, int size) {
    return (fabs(sum) - size * meanScore) / pow(size, scalingExponent);
  }
  
 public:
  // read graph from file
  void read(int nn, const string& edgeFile, map<string,int>&,
	    const vector<string>&);

  // structural queries
  // isEdge is slow, but only currently used for output, so thats ok
  int degree(int node) { return _degree[node]; }
  vector<int> ball(int center, int radius); // returns B(center, radius)

  // main algorithm
  void setScores(vector<double>&, double, double);
  double findSubgraph(int root, const Args& arg, vector<int>& best);
  // output
  void print(void); // used mostly for debugging
  ostream& writeDOT(ostream&, const SubGraph&); // write graph as DOT
};

#endif
