#include "Graph.h"
#include "Data.h"

inline bool vectorFind(const int key, const vector<int>& vec) {
  for (int i = 0; i < vec.size(); i++)
    if (vec[i] == key)
      return true;
  return false;
}

void Graph::addEdge(int v1, int v2, string type, string source) {
  pair<int,int> edge(v1,v2);
  if (edgeType.find(edge) == edgeType.end()) { // new edge
    if (!vectorFind(v2, neighbors[v1])) // wasteful search but fast enough
      neighbors[v1].push_back(v2);
    if (!vectorFind(v1, neighbors[v2]))
      neighbors[v2].push_back(v1);
    edgeType[edge] = type;
    edgeSource[edge] = source;
    nEdges++;
  }
  else { // do nothing for now, may want to track the various types and sources
  }
}

// correct but wasteful; should declare/reset newBoundary inside the for loop
vector<int> Graph::ball(int center, int radius) {
  vector<int> taken(nNodes); // tracks nodes already in, avoids binary search
  taken[center] = 1;
  vector<int> ball, boundary, newBoundary;
  ball.push_back(center);
  boundary.push_back(center);

  vector<int>::iterator vi;
  vector<int>::iterator wi;
  int depth;
  for (depth = 1; depth <= radius; depth++) { // build ball one layer at a time
    for (wi = boundary.begin(); wi != boundary.end(); wi++)
      for (vi = neighbors[*wi].begin(); vi != neighbors[*wi].end(); vi++)
	if (taken[*vi] == 0) {
	  ball.push_back(*vi);
	  newBoundary.push_back(*vi);
	  taken[*vi] = 1;
	}
    boundary = newBoundary;
  }
  return ball;
}

// the surrogate constructor
//
// each line in the edgeFile gives the two endpoints of an edge, followed
// optionally by an interaction type and source
//
// if either endpoint is not found in id2gene, line is silently ignored
// loops (vertex connected to itself) are ignored
void Graph::read(int nn, const string& edgeFile, map<string,int>& id2gene,
		 const vector<string>& label) {
  nNodes = nn; // set basic variables
  nEdges = 0;
  _degree.resize(nNodes);
  neighbors.resize(nNodes);

  if (label.size() > 0) // set labels
    _label = label; // should make sure label.size() == nNodes !!
  else {
    _label.resize(nNodes);
    for (int i = 0; i < nNodes; i++)
      _label[i] = i; // default label is the integer value of the node
  }

  ifstream fs(edgeFile.c_str());
  string str;
  while (getline(fs, str)) {
    istringstream ss(str);
    string v1, v2, itype = "", isource = "";
    ss >> v1 >> v2 >> itype >> isource;
    map<string,int>::iterator m1 = id2gene.find(v1), m2 = id2gene.find(v2);
    if (v1 != v2 && m1 != id2gene.end() && m2 != id2gene.end()) {
      addEdge(m1->second, m2->second, itype, isource);
    }
  }
  for (int i = 0; i < nNodes; i++) // set degrees
    _degree[i] = neighbors[i].size();
  cerr << "Read " << nEdges << " edges from " << edgeFile << endl;
}

class CompareScores {
  const vector<double> *_pscore;
 public:
  CompareScores(vector<double> *pscore) : _pscore(pscore) {}
  bool operator()(const int x, const int y) const {
    return (*_pscore)[x] < (*_pscore)[y];
  }
};

void Graph::setScores(vector<double>& T, double mean, double sd) {
  // static int nCalls = 0; // used for debugging
  score = T;
  meanScore = mean;
  sdScore = sd;
  CompareScores comp(&score);
  for (int i = 0; i < nNodes; i++) // sort neighbors in order of score
    sort(neighbors[i].begin(), neighbors[i].end(), comp);
}

// main algorithm; must call setScores before using it
// finds a high-scoring subgraph, returns its score
// NEEDS TESTING if scores can be negative; MAY FAIL
double Graph::findSubgraph(int root, const Args& arg, vector<int>& best) {
  int depth = arg.depth;
  scalingExponent = arg.scalingExponent;
  int i, j, k, m;
  vector<int>::iterator vi; 
  static vector<bool> taken(nNodes); // marks nodes in current subgraph; static gives slight speedup
  double currentSum, currentScore;
  double maxScore; // best subgraph score found so far
  vector<int> subgraph; // best subgraph found so far
  vector<int> subdegree; // degrees for subgraph spanning tree, used by Metropolis step
  vector<int> parent; // index of parent in spanning tree
  
  // build subgraph by greedy expansion
  for (i = 0; i < 1; i++) { // superfluous loop, may be useful in future randomized version
    int node, bestNode, bestk;
    vector<int> index1, index2;
    subgraph.reserve(depth); // in future randomized version, must reset subgraph
    index1.reserve(depth);
    index2.reserve(depth);
    subdegree.reserve(depth);
    parent.reserve(depth);
    subgraph.push_back(root); // add root to relevant arrays
    index1.push_back(0);
    index2.push_back(neighbors[root].size() - 1);
    subdegree.push_back(0);
    parent.push_back(-1); // root has no parent
    taken[root] = true;
    currentSum = score[root]; // sum of scores for current subgraph
    currentScore = sum2score(currentSum, 1); // currentSum, scaled
    // increase subgraph, one vertex at a time
    for (j = 1; j < depth; j++) {
      double bestSum = currentSum;
      bestNode = -1;
      if (true) { // pick node with max score; use sums to handle scores < 0
        for (k = 0; k < subgraph.size(); k++) { // for each vertex, find the best neighbor
          node = subgraph[k];
          for (m = index2[k]; m >= 0 && taken[neighbors[node][m]]; m--); // skip taken vertices
          index2[k] = m;
          if (m >= 0) {
            double s = currentSum + score[neighbors[node][m]];
            if (fabs(s) > fabs(bestSum)) { // this node is better
              bestSum = s;
              bestNode = neighbors[node][m];
              bestk = k;
            } // max so far found
          }
        }
      }
      else { // randomized search, not implemented
      }
      if (bestNode == -1) // there are no nodes to be added
        break;
      double newScore = sum2score(bestSum, j+1);
      if (arg.flexSize && newScore < currentScore) // node not good enough
        break; // finish the search, do not add node
      else {
        subgraph.push_back(bestNode); // add node
        index1.push_back(0);
        index2.push_back(neighbors[bestNode].size() - 1);
        subdegree.push_back(1);
        subdegree[bestk]++;
        parent.push_back(bestk);
        taken[bestNode] = true;
        currentSum = bestSum;
        currentScore = newScore;
      }
    } // finished greedy expansion
    if (i == 0 || maxScore < currentScore) {
      best = subgraph;
      maxScore = currentScore;
    }
  } // finished superfluous loop
  
  // now try Metropolis updates
  // if changing algorithm, make sure currentSum is correct
  bool DEBUG_MET = false;
  double metStartScore = currentScore; // used for debugging
  double kMet = arg.kMet / sdScore; // small makes it easy to jump
  int tries = 0, moves = 0, maxes = 0;
  for (k = 0; k < arg.iterations; k++) {
    int v1 = urand(subgraph.size() - 1) + 1; // root must remain in subgraph
    int v2 = urand(subgraph.size());
    int node1 = subgraph[v1]; // random nodes in subgraph
    int node2 = subgraph[v2];
    if (subdegree[v1] == 1 && v1 != v2) {
      int v3 = urand(neighbors[node2].size());
      for (vi = neighbors[node2].end() - 1; vi != neighbors[node2].begin() && taken[*vi]; vi--);
        // fix: should choose random neighbor of node2
      if (!taken[*vi] && node1 != *vi) {
        int node3 = *vi;
        ++tries;
        double scoreDiff = score[node3] - score[node1];
        double threshold = exp(scoreDiff * kMet);
        if (DEBUG_MET) {
          cout.precision(3);
          cout << root << ":" << k << "[" << currentScore << "] ";
          cout << node1 << "[" << score[node1] << "," << subdegree[v1] << "] ";
          cout << node2 << "[" << score[node2] << "," << subdegree[v2] << "] ";
          cout << node3 << "[" << score[node3] << "] ";
          cout << threshold << " ";
        }
        if (threshold > rand() / (double) RAND_MAX) { // perform update
          ++moves;
          currentSum += scoreDiff;
          currentScore = sum2score(currentSum, j);
          taken[node3] = true;
          taken[node1] = false;
          subgraph[v1] = node3;
          subdegree[v2]++; // node2 now has one extra child
          subdegree[parent[v1]]--; // parent of node1 lost a child
          parent[v1] = v2;
          if (currentScore > maxScore) {
            ++maxes;
            best = subgraph;
            maxScore = currentScore;
          }
        }
        if (DEBUG_MET)
          cout << moves << "/" << tries << endl;
      }
    }
  } // done with Metropolis
  if (DEBUG_MET) {
    cout << "METROPOLIS " << root << ' ' << maxes << "/" << moves << "/" << tries << " ";
    cout << metStartScore << "/" << maxScore << endl;
  }
  // since taken is static, must restore it to all false
  for (k = 0; k < subgraph.size(); k++)
    taken[subgraph[k]] = false;
  return maxScore;
}

void Graph::print(void) {
  int i;
  vector<int>::iterator si;
  cout << nNodes << " nodes " << nEdges << " edges\n";
  for (i = 0; i < nNodes; i++) // for each node
    if (_degree[i] != 0) {
      cout << i << ' ' << _degree[i] << " :"; // print node, degree
      for (si = neighbors[i].begin(); si != neighbors[i].end(); si++)
	      cout << ' ' << *si; // then print all neighbors
      cout << endl;
    }
}

// helper functions for writeDOT
// note that v1, v2 are the internal ids, NOT the Entrez ids
string getDOTEdge(const int v1, const int v2, const string type) {
  stringstream edge;
  if (type == "I") // undirected interaction
    edge << v1 << " -> " << v2 << " [arrowhead = \"none\" ]";
  else
    edge << v1 << " -> " << v2 << " [label = \"" << type << "\" ]";
  return edge.str();
}

string getDOTNode(const int v, const string label, const string color = "") {
  stringstream node;
  string cstring = (color == "" ? "" : ", color=\"" + color + "\"");
  node << v << " [label = \"" << label << "\"" << cstring << " ] ; \n";
  return node.str();
}

ostream& Graph::writeDOT(ostream& str, const SubGraph& sg) {
  bool PrintContext = false; // if true, prints all nearest neighbors
  int i, j, node;
  str << "digraph G {\n"; // print header
  str << "overlap = scale ;\n";
  
  set<int> printedNodes; // keep track of nodes printed to avoid duplication
  for (i = 0; i < sg.size(); i++) // print nodes in sg
    if (printedNodes.find(sg[i]) == printedNodes.end()) { // not yet printed
      str << getDOTNode(sg[i], _label[sg[i]]);
      printedNodes.insert(sg[i]);
    }
  if (PrintContext) // print neighbors to sg
    for (i = 0; i < sg.size(); i++)
      for (j = 0; j < neighbors[sg[i]].size(); j++) {
        node = neighbors[sg[i]][j];
        if (printedNodes.find(node) == printedNodes.end()) { // not yet printed
          str << getDOTNode(node, _label[node]);
          printedNodes.insert(node);
        }
      }
  
  set< pair<int,int> > printEdges; // edges to print
  set< pair<int,int> >::iterator si;
  for (i = 0; i < sg.size(); i++)
    for (j = 0; j < neighbors[sg[i]].size(); j++) {
      pair<int,int> edge1(sg[i], neighbors[sg[i]][j]);
      pair<int,int> edge2(neighbors[sg[i]][j], sg[i]);
      if (PrintContext) { // print all edges from i
        printEdges.insert(edge1);
        printEdges.insert(edge2);
      }
      else if (vectorFind(neighbors[sg[i]][j], sg)) // print only edges in sg
        printEdges.insert(edge1); // edge2 will be inserted when processing j
    }
  // now print edges; each edge appears twice in printEdges, so we need to check
  // edgeType to see which direction to print (possibly both)
  for (si = printEdges.begin(); si != printEdges.end(); si++)
    if (edgeType.find(*si) != edgeType.end())
      str << getDOTEdge(si->first, si->second, edgeType[*si]) << endl;

  str << "}\n"; // print trailer
  return str;
}
