// class for gene expression data

#include "GEData.h"
#include "cdf.h"

// input functions
// the constructor first reads the probe2entrez map; it also assigns an
//   integer id to each gene based on order in the file
// then it reads gene names and phenotypes
// then it reads expression
//   probes that dont map to a gene (in probe2entrez) are ignored
//   genes with no probes get expression 0
//   genes with multiple probes get probe average
// then it reads the interaction graph
//   genes not in probe2entrez get ignored

int GEData::readMap(const char *file) {
  ifstream fs(file);
  int lines(0), n(0);
  string str;
  while (getline(fs, str)) {
    lines++;
    istringstream ss(str);
    string probe, id = "NA";
    ss >> probe >> id;
    if (ss && probe != "" && id != "NA") { // good line
      if (probe2id.find(probe) == probe2id.end()) { // new probe
	probe2id[probe] = id;
	if (id2gene.find(id) == id2gene.end()) { // new gene
	  id2gene[id] = n;
	  geneID.push_back(id);
	  n++;
	}
      }
      else { // if probe appeared before, check it has same id
	if (probe2id[probe] != id) {
	  cerr << "readMap: multiple IDs ";
	  cerr << probe << ' ' << probe2id[probe] << ' ' << id << endl;
	}
      }
    }
  }
  if (n == 0) {
    cerr << "Bad or missing map file " << file << endl;
    exit(1);
  }
  cerr << "readMap: read " << lines << " lines and ";
  cerr << n << " genes from " << file << endl;
  return n;
}

int GEData::readPhenotype(const char *file) {
  ifstream fs2(file);
  string name;
  nPheno = 0;
  while (fs2 >> name) {
    if (phenoMap.find(name) == phenoMap.end()) { // new type
      n2pheno.push_back(name);
      phenoMap[name] = nPheno++;
    }
    origPhenotype.push_back(phenoMap[name]);
  }

  int k = origPhenotype.size();
  phenotype.resize(k);
  if (k == 0) {
    cerr << "Bad phenotype file " << file << ". Quitting.\n";
    exit(1);
  }
  else if (k == 1) {
    cerr << "Found only one phenotype, need at least two.\n";
    exit(1);
  }
  else {
    cerr << "Read " << nPheno << " phenotypes for " << k
	 << " arrays from " << file << ": ";
    for (int i = 0; i < nPheno; i++)
      cerr << n2pheno[i] << ' ';
    cerr << endl;
  }

  // now that we know the phenotypes, decide what the score should be
  if (arg.scoreType == SCORE_UNKNOWN) {
    if (nPheno > 2 && isAll(arg.compare)) { // use F statistic
      if (arg.algoType == ALGO_GXNA || arg.radius == 0) {
	arg.scoreType = SCORE_SIGMA_F;
	cerr << "Using score Sigma(F)\n";
      }
      else {
	arg.scoreType = SCORE_F_SIGMA;
	cerr << "Using score F(Sigma)\n";
      }
    }
    else { // use T statistic
      if (arg.algoType == ALGO_GXNA || arg.radius == 0) {
	arg.scoreType = SCORE_SIGMA_T;
	cerr << "Using score Sigma(T)\n";
      }
      else {
	arg.scoreType = SCORE_T_SIGMA;
	cerr << "Using score T(Sigma)\n";
      }
    }
  }

  if (isT(arg.scoreType) && k >= 2) {
    if (isAll(arg.compare))
      cerr << "No phenotypes to compare specified, will compare "
	   << n2pheno[0] << " to " << n2pheno[1] << endl;
    else {
      int i = arg.compare.find(",");
      if (i == string::npos) {
        cerr << "Bad argument -compare: " << arg.compare << endl;
	exit(1);
      }
      string p1 = arg.compare.substr(0, i);
      string p2 = arg.compare.substr(i+1);
      if (phenoMap.find(p1) == phenoMap.end()) {
	cerr << "Phenotype not found: " << p1 << endl;
	exit(1);
      }
      if (phenoMap.find(p2) == phenoMap.end()) {
	cerr << "Phenotype not found: " << p2 << endl;
	exit(1);
      }
      if (p1 == p2) {
	cerr << "Cannot compare phenotype " << p1 << " to itself\n";
	exit(1);
      }
      FastDataSet::setPheno(phenoMap[p1], phenoMap[p2]);
    }
  }
  return k;
}

// when multiple probes on an array map to the same gene, readExpression()
// averages their expression values [there may be better ways]
//
// need to improve error checking. currently:
//   if a line is too long then the extra values are ignored
//   missing values are not allowed
int GEData::readExpression(const char *file) {
  ifstream fs(file);
  int good(0), bad(0), i;
  string str;
  double mm;
  vector<FastDataSet> mvec(nGenes * nArrays);
  probesPerGene.resize(nGenes);
  while (getline(fs, str)) {
    string probe;
    istringstream ss(str);
    ss >> probe;
    map<string,string>::iterator pp = probe2id.find(probe);
    if (ss && pp != probe2id.end()) { // probe maps to known gene
      int gene = id2gene[pp->second];
      ++probesPerGene[gene];
      int index = gene * nArrays;
      for (i = 0; i < nArrays; i++) {
	if (ss >> mm)
	  mvec[index+i].insert(mm);
	else {
	  cerr << "Bad expression line " << str << "\nQuitting.\n";
	  exit(1);
	}
      }
      good++;
    }
    else bad++; // dont know how to map probe to a gene
  }

  for (i = 0; i < nGenes * nArrays; i++)
    M[i] = mvec[i].mean();

  if (good + bad == 0) {
    cerr << "Bad expression file " << file << ". Quitting.\n";
    exit(1);
  }
  else
    cerr << "readExpression: " << good << " good lines, " 
	 << bad << " skipped lines in " << file << endl;
  return good;
}

// constructor

GEData::GEData(Args& _arg) :
  arg(_arg),
  nTypes(0),
  permCount(0),
  meanT(0),
  gbest(arg.depth),
  DOTPrinted(0),
  baseFileName(arg.name + "_" + arg.version),
  outputFileName(arg.outputDir + "/" + baseFileName + ".res"),
  htmlFileName(arg.outputDir + "/" + baseFileName + ".html"),
  htmlFrameFileName(arg.outputDir + "/" + baseFileName + "_frame1.html"),
  frameStream(htmlFrameFileName.c_str())
{
  // read id2name map
  map<string,string> id2name;
  string id, name;
  int i, k = 0;
  ifstream fs1(arg.namesFile.c_str());
  while (fs1 >> id >> name)
    id2name[id] = name;

  // read probe2id map and phenotype
  nGenes = readMap(arg.mapFile.c_str());
  nArrays = readPhenotype(arg.phenoFile.c_str());
  geneName.resize(nGenes);
  for (i = 0; i < nGenes; i++) {
    if (id2name.find(geneID[i]) != id2name.end())
      geneName[i] = id2name[geneID[i]];
    else
      geneName[i] = "NA";
  }
  groupScore = new double[nArrays];
  M = new double[nGenes * nArrays];

  // read secondary phenotype, if needed
  if (arg.useInvariantPerms) {
    ifstream fs3(arg.typesFile.c_str());
    map<string,int> tmap; // maps type names to integers
    arrayType.resize(nArrays);
    while (fs3 >> name && k < nArrays) {
      if (tmap.find(name) == tmap.end()) // new type
	tmap[name] = nTypes++;
      arrayType[k] = tmap[name];    
      ++k;
    }
    if (k < nArrays) {
      cerr << "Bad types file " << arg.typesFile << ". Quitting.\n";
      exit(1);
    }
    else
      cerr << "Read " << k << " arrays, " << nTypes << " cell types from " << 
	arg.typesFile << ' ' << endl;
  }

  readExpression(arg.expFile.c_str());

  // set some more parameters
  if (arg.algoType == ALGO_GXNA) {
    if (isT(arg.scoreType))
      arg.takeAbs = true; // otherwise graph search may fail
    if (arg.scoreType == SCORE_T_SIGMA) // graph search needs sigmaT or sigmaF
      arg.scoreType = SCORE_SIGMA_T;
    else if (arg.scoreType == SCORE_F_SIGMA)
      arg.scoreType = SCORE_SIGMA_F;
  }
  else if (arg.radius == 0) // single gene analysis
    arg.takeAbs = false;

  T.resize(nGenes); // T will be computed by newPermutation
  signT.resize(nGenes);
  sd.resize(nGenes);
  foldChange.resize(nGenes);
  sfac.resize(nGenes);
  FastDataSet lsd; // tracks log of gene sds (used for shrinkage)
  FastDataSet gvar; // gene variances
  double dfData = nArrays - nPheno; // degrees of freedom in data

  for (i = 0; i < nGenes; i++) { // compute global SDs, fold change etc
    FastDataSet d, d0, d1;
    for (int j = 0; j < nArrays; j++) {
      double x = M[i * nArrays + j];
      d.insert(x);      
      if (origPhenotype[j] == FastDataSet::pheno0)
	d0.insert(x);
      else if (origPhenotype[j] == FastDataSet::pheno1)
	d1.insert(x);
    }
    foldChange[i] = d1.mean() - d0.mean();
    sd[i] = d.sd();
    if (sd[i] > 0) { // constant probes are presumably bad so we ignore them
      lsd.insert(2 * log(sd[i]));
      gvar.insert(sd[i] * sd[i]);
    }
  }

  // estimate parameters for shrinkage factors; sloppy but simple
  double shrinkDF = arg.shrinkDF;
  double shrinkVar = arg.shrinkVar;
  if (arg.shrink && nGenes > 1 && dfData > 0) { // re-estimate the parameters
    double yy = lsd.var() * (1 + 1 / (nGenes - 1.0)) - trigamma(dfData / 2.0);
    double d0, var0;
    if (yy > 0) {
      d0 = 2 * trigammainv(yy);
      var0 = exp(lsd.mean() + digamma(d0 / 2) + log(dfData / d0) - digamma(dfData / 2));
    }
    else {
      d0 = 100 * dfData; // just need something large
      var0 = exp(lsd.mean() + log(dfData / 2) - digamma(dfData / 2));
    }
    shrinkDF = d0;
    shrinkVar = var0;
  }
  cerr << "Gene variances: mean = " << gvar.mean() << " sd = " << gvar.sd() << endl;
  if ((arg.scoreType == SCORE_T_SIGMA || arg.scoreType == SCORE_F_SIGMA) 
      && arg.shrinkDF > 0) {
    arg.shrinkDF = 0;
    cerr << "No shrinkage for this scoreType, setting shrinkDF to 0\n";
  }
  if (shrinkDF > 0)
    cerr << "Shrinkage parameters: df = " << shrinkDF << " var = " << shrinkVar << endl;

  for (i = 0; i < nGenes; i++) { // compute shrinkage factors
    double sf = 1; // by default, no shrinkage
    if (shrinkDF > 0 && sd[i] > 0) { // compute shrinkage factors
      if (dfData > 0) { // avoid division by zero
	double dfRatio = shrinkDF / dfData;
	sf = (dfRatio + 1) / (dfRatio * shrinkVar / sd[i] / sd[i] + 1);
      }
    }
    sfac[i] = sqrt(sf);
  }

  vector<string> label(nGenes); // labels for graph nodes
  for (i = 0; i < nGenes; i++) { // label contains name and Tstat
    ostringstream slab;
    slab.precision(2);
    double x = (isT(arg.scoreType) ? foldChange[i] : origScore(i));
    slab << geneName[i] << "\\n" << int(x * 100) / 100.0;
    label[i] = slab.str();
  }
  gg.read(nGenes, arg.edgeFile, id2gene, label); // read the interaction graph
  subGraphs.resize(nGenes);
  subGraphScore.resize(nGenes);
}

// master computation functions
// testBall and testSubgraphs share the same output function
// in order to enable this, there is some redundancy in testBall
// also, this relies on SubGraph being vector<int> - should fix...

void GEData::testBall(PermutationGenerator& pg) {
  vector<int> vec; // stores the roots to be tested
  for (int gene = 0; gene < nGenes; gene++) // filter the roots
    if (gg.degree(gene) >= arg.minDegree && sd[gene] >= arg.minSD) {
      subGraphs[gene] = gg.ball(gene, arg.radius);
      vec.push_back(gene);
    }
  MultipleTest<int,GEData> mt(*this, vec);
  if (arg.maxTscaled) // multiple test
    mt.maxTscaled(pg);
  else
    mt.maxT(pg);
  printResults(mt);
}

// potential roots will not be within arg.radius of each other
void GEData::testSubgraphs(PermutationGenerator& pg) {
  vector<int> vec; // list of roots for subgraph search
  vector<bool> taken(nGenes); // false if gene is a potential root
  for (int gene = 0; gene < nGenes; gene++) // filter the roots
    if (gg.degree(gene) >= arg.minDegree && sd[gene] >= arg.minSD && 
	!taken[gene]) {
      vec.push_back(gene); // add gene as a root
      vector<int> ball = gg.ball(gene, arg.radius);
      for (int i = 0; i < ball.size(); i++)
	taken[ball[i]] = true; // dont use nearby genes as roots
    }
  MultipleTest<int,GEData> mt(*this, vec);
  if (arg.maxTscaled) // multiple test
    mt.maxTscaled(pg);
  else
    mt.maxT(pg);
  printResults(mt);
}

// groupSum depends on current perm only if takeAbs is true
double GEData::groupSum(const vector<int>& genes, int array) {
  double sum = 0;
  vector<int>::const_iterator si;
  for (si = genes.begin(); si != genes.end(); si++) {
    double m = M[(*si) * nArrays + array];
    if (arg.takeAbs && signT[*si] < 0)
      sum -= m;
    else
      sum += m;
  }
  return sum;
}

// main low level scoring function
double GEData::vectorScore(double *values, const vector<int>& pheno) {
  if (isF(arg.scoreType)) { // F stat
    double f = fstat(values, pheno, nPheno);
    return zfCDF(f, nPheno - 1, pheno.size() - nPheno);
  }
  else { // T stat; should convert to z-score
    return tstatPheno(values, pheno);
  }
}

// housekeeping for a new permutation
void GEData::newPermutation(int *perm) {
  for (int i = 0; i < nArrays; i++)
    phenotype[i] = origPhenotype[perm[i]]; // compute new phenotype
  meanT = 0;
  double t;
  FastDataSet tset;
  for (int j = 0; j < nGenes; j++) { // recompute T and meanT
    t = vectorScore(M + j * nArrays, phenotype);
    if (arg.shrinkDF > 0)
      t *= sfac[j];
    T[j] = (arg.takeAbs ? fabs(t) : t);
    signT[j] = (t < 0 ? -1 : 1); // ignore case t == 0, should be ok
    tset.insert(T[j]);
  }
  meanT = (arg.takeAbs ? tset.mean() : 0); // no need to center if Tstats are signed
  gg.setScores(T, meanT, tset.sd());
  ++permCount;
}

// scoring function for MultipleTest, predefined case
double GEData::scoreSubGraph(const SubGraph& genes) {
  if (arg.scoreType == SCORE_SIGMA_T || arg.scoreType == SCORE_SIGMA_F) {
    double sum = 0;
    int s = genes.size();
    for (int i = 0; i < s; i++) // use values stored by newPermutation()
      sum += T[genes[i]];
    return (sum - s * meanT) / pow(s, arg.scalingExponent);
  }
  else { // T(sum) or F(sum)
    for (int i = 0; i < nArrays; i++) // compute group scores
      groupScore[i] = groupSum(genes, i);
    return vectorScore(groupScore, phenotype);
  }
}

double GEData::operator() (int n) {
  double score;
  if (arg.algoType == ALGO_GXNA) { // adapted case
    // find best graph rooted at gene n, store it in gbest
    score = gg.findSubgraph(n, arg, gbest); // findSubgraph guarantees score >= 0
    if (permCount == 1) { // for the real data (id perm), save the best graph
      subGraphs[n] = gbest;
      subGraphScore[n] = score;
    }
    return score;
  }
  else { // predefined group
    score = scoreSubGraph(subGraphs[n]);
    if (permCount == 1)
      subGraphScore[n] = score;
  }
  return fabs(score);
}

// output functions

// type should be "dot" or "txt" or "svg"
string GEData::graphFileName(int i, const string& type) {
  ostringstream fileName;
  fileName << baseFileName << "_" << i << "." << type; 
  return fileName.str();
}

// main output function
void GEData::printResults(const MultipleTest<int,GEData>& mt) {
  beginHTML(); // prepare html files
  vector<bool> printed(nGenes); // initialized to all false
  int i, j, rows = 0;
  ofstream fs(outputFileName.c_str());
  fs.precision(4);
  for (i = 0; i < mt.size(); i++) {
    int gene = mt.getObject(i); // print to output file
    fs << i << ' ';
    print(fs, gene); // print gene
    double rawP = mt.getRawP(i), adjP = mt.getAdjP(i);
    fs << ' ' << rawP << ' ' << adjP << endl;

    // now check if we print to html and make graphs
    int n1 = 0, n2 = subGraphs[gene].size();
    for (j = 0; j < n2; j++)
      if (printed[subGraphs[gene][j]])
	n1++;
    if (n1 <= arg.maxOverlap * n2) { // ok to print
      for (j = 0; j < n2; j++)
	printed[subGraphs[gene][j]] = true; // mark as printed
      if (DOTPrinted < arg.maxRows)
	addRow(i, geneName[gene], geneID[gene], subGraphs[gene].size(),
	       subGraphScore[gene], rawP, adjP);
      if (DOTPrinted < arg.graphCount)
	makeDOT(gene, DOTPrinted);
      DOTPrinted++;
    }
  }
  endHTML();
}

ostream& GEData::print(ostream& str, int gene) {
  str << geneInfo(gene) << "  "; // print root details
  for (int i = 0; i < subGraphs[gene].size(); i++)
    str << geneID[subGraphs[gene][i]] << ' '; // print IDs of nodes in graph
  str << subGraphScore[gene];
  return str;
}

// make DOT, TXT, possibly SVG
void GEData::makeDOT(int gene, int nDOT) {
  string DOTFile = arg.outputDir + "/" + graphFileName(nDOT, "dot");
  ofstream sdot(DOTFile.c_str());
  gg.writeDOT(sdot, subGraphs[gene]); // write to file
  sdot.close();

  string TXTFile = arg.outputDir + "/" + graphFileName(nDOT, "txt"); // make TXT file
  ofstream stxt(TXTFile.c_str());
  stxt.precision(3);
  for (int i = 0; i < subGraphs[gene].size(); i++) {
    int v = subGraphs[gene][i];
    stxt << setw(10) << geneName[v] << ' ';
    stxt << setw(6) << geneID[v] << ' ';
    stxt << setw(2) << probesPerGene[v] << ' ';
    stxt << setw(7) << foldChange[v] << ' ';
    stxt << setw(7) << sd[v] << ' ';
    stxt << setw(7) << origScore(v) << endl;
  }
  stxt.close();

  if (arg.runDOT) { // make SVG
    string SVGFile = arg.outputDir + "/" + graphFileName(nDOT, "svg");
    string cmd = "neato -Tsvg " + DOTFile + " >" + SVGFile;
    system(cmd.c_str());
    /* SVG change currently off
    // now we need to change the <SVG header line in the SVG file
    // this will make it display better in web browsers
    ifstream svg1(SVGFile.c_str()); // read into memory
    vector<string> svgvec;
    string str;
    while (getline(svg1, str))
      svgvec.push_back(str.substr(0, 4) == "<svg" ? "<svg" : str);
    svg1.close();
    ofstream svg2(SVGFile.c_str()); // write back to disk
    for (int i = 0; i < svgvec.size(); i++)
      svg2 << svgvec[i] << endl;
    svg2.close();
    */
  }
}

void GEData::beginHTML() {
  string startingFrame;
  if (arg.runDOT) // starting frame is SVG file
    startingFrame = baseFileName + "_0.svg";
  else
    startingFrame = baseFileName + "_0.txt";

  // write main html file
  ofstream html(htmlFileName.c_str());
  html << "<html>" << endl;
  html << "<title>" << "GXNA " << arg.name << ' ' << arg.version << "</title>" << endl;
  html << "<frameset cols=\"35%,65%\">" << endl;
  html << "<frame src=\"" << baseFileName << "_frame1.html" << "\">" << endl;
  html << "<frame src=\"" << startingFrame << "\" name=\"frame2\">" << endl;
  html << "</frameset>" << endl;
  html << "</html>" << endl;

  // start frame html file
  frameStream << "<html>" << endl;
  frameStream << "<table border=\"1\">" << endl;
  frameStream << "<tr>" << endl;
  frameStream << "<th>rank</th>" << endl;
  frameStream << "<th>root</th>" << endl;
  frameStream << "<th>rootid</th>" << endl;
  frameStream << "<th>size</th>" << endl;
  frameStream << "<th>score</th>" << endl;
  frameStream << "<th>rawp</th>" << endl;
  frameStream << "<th>adjp</th>" << endl;
  frameStream << "</tr>" << endl;
  frameStream.precision(3);
}

void GEData::addRow(int n, const string& root, const string& rootid, int size,
		    double score, double rawp, double adjp) {
  frameStream << "<tr>" << endl;
  if (DOTPrinted < arg.graphCount) {
    string url;
    if (arg.runDOT)
      url = graphFileName(DOTPrinted, "svg");
    else
      url = graphFileName(DOTPrinted, "txt");
    frameStream << "<td><a href=\"" << url << "\" target=\"frame2\">" << n << "</a></td>" << endl;
  }
  else
    frameStream << "<td>" << n << "</td>" << endl;
  frameStream << "<td>" << root << "</td>" << endl;
  frameStream << "<td>" << rootid << "</td>" << endl;
  frameStream << "<td>" << size << "</td>" << endl;
  frameStream << "<td>" << score << "</td>" << endl;
  frameStream << "<td>" << rawp << "</td>" << endl;
  frameStream << "<td>" << adjp << "</td>" << endl;
  frameStream << "</tr>" << endl;
}

void GEData::endHTML() {
  frameStream << "</table>" << endl;
  frameStream << "</html>" << endl;
  frameStream.close();
}
