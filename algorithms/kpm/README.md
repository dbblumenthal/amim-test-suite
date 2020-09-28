
Standalone version of KeyPathwayMiner 5.0
============
Given a biological network and a set of case-control studies, KeyPathwayMiner efficiently extracts all maximal connected sub-networks. These sub-networks contain the genes that are mainly dysregulated, e.g., differentially expressed, in most cases studied.

For more information please visit our website [(Key Pathway Miner website)](https://exbio.wzw.tum.de/keypathwayminer/).

About
=================
    KeyPathwayMiner version 5.0
    Copyright 2016 by
    Nicolas Alcaraz: nalcaraz@binf.ku.dk 
    Jan Baumbach: jan.baumbach@wzw.tum.de
    Markus List: markus.list@wzw.tum.de
    Standalone version adapted for KPM 5.0: Konstantinos Mechteridis

Overview
=================
<!--ts-->
   * [About](#about)
   * [Usage](#usage)
   * [Input files format](#input-files-format)
   * [Parameters](#parameters)
      * [Input file parameters](#input-file-parameters)
      * [Input options basic parameters](#input-options-basic-parameters)
      * [Output parameters](#output-parameters)
      * [Advanced options](#advanced-options)
<!--te-->

Usage
=================
   General structure:
      
    java -jar [jvm options] KPM-5.0.jar [-KEY1=VAL1] .... [-KEYN=VALN] 

   Simple executions examples:
   
    java -jar -Xmx2G KPM-5.0.jar -strategy=INES -algo=GREEDY -K=2 -L1=5
      
    java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5
      
    java -jar KPM-5.0.jar -strategy=GLONE -algo=ACO -datasetsFile=resources/datasets_file.txt

   Combine multiple datasets:
   
    java -jar KPM-5.0.jar -strategy=GLONE -algo=GREEDY -L1=5 -matrix1=resources/datasets/colon-gene-expression-DOWN-p0.05.txt -L2=6 -matrix2=resources/datasets/colon-gene-expression-UP-p0.05.txt
    
   Use ranged values with batch:
   
    java -jar KPM-5.0.jar -batch  -L1_batch=1,2,3 -L2_batch=1,2,3 -K_batch=1,2,3 -strategy=INES -algo=GREEDY -matrix1=resources/datasets/colon-gene-expression-DOWN-p0.05.txt -matrix2=resources/datasets/colon-gene-expression-UP-p0.05.txt 
    
   Use perturbation:
   
    java -jar KPM-5.0.jar -strategy=INES -algo=GREEDY -K=2 -L1=5 -perturbation=10,10,20,1 -perturbation_technique=nodeswap
        
   For Help:

    java -jar KPM-5.0.jar -help

   Note: If the input is large and/or complex then the java virtual machine options must be set.

Input files format
=====================

    KPM takes as input several files

    -- GRAPH file:  The file containing all the interactions
    of the protein interaction network.  This must be in sif format:

    NODE1   INTERACTION_TYPE   NODE2
    112 pp 342
    12  pp 342
    ...
    ...

    -- MATRIX file(s):  The file containing and indicator matrix
    for expression studies. KPM can take as input several of these
    files, which can be defined either through the command line 
    or in the DATASETS file (see below). 

    GENE_ID   CASE1	  CASE2   CASE3
    10203         1	      0       1
    3232          0           0       1
    ...
    ...		

    -- DATASETS file:  This file contains the paths to each individual
    indicator matrix file and it's corresponding L parameter. The format
    should be the following:

    ID      L       PATH
    1       10      path/to/matrix1.txt
    2       15      path/to/matrix2.txt
    ...
    ...

    -- POSITIVE/NEGATIVE list file: these are optional files that contain
    a list of genes that will be given high (POSITIVE list) or low (NEGATIVE list)
    priority when searching for pathways.  The format is:

    NODE_ID
    11423
    1213
    ...
    ...

Parameters
=================

    All parameters can be defined in the "kpm.properties" file
    and most through command line arguments. In case the same parameter
    is defined both in the properties file and through a command
    line argument, then the command line argument will have preference
    and override the value in the properties file. 

Input file parameters
-----
### Graph File

      -graphFile {string}   
       The path to the graph file. Must be in sif format

      -gfHeader {boolean}
       If the graph file has a header row 

      -gfSep {string}
       The separating character for the columns
       in the graph file: TAB, SPACE or COMMA 

     -isDirected {boolean}
       If graph should be considered directed or not

### Matrix File

      -matrixN {string}
       The path to the n-th matrix 

      -mfHeader {boolean}
       If the matrix file has a header row 

      -mfSep {String}
       The separating character for the columns
       in the matrix file: TAB, SPACE or COMMA 

### Datasets File

      -datasetsFile {string}  
       Path to the file containing the list of
       paths to the matrices and their respective case exceptions
       (L parameters) 

      -dfHeader {boolean}
       If the datasets file contains a header

      -dfSep {string}
       The separating character for the columns
       in the datasets file: TAB, SPACE or COMMA 

### Positive & negative list file

      -positiveFile  {string}
       The path to the file with the positive gene list
  
      -negativeFile {string} 
        The path to the file with the negative gene list

### Validation file
     -validationFile {string}
      Gold-standard set to determine how enriched extracted pathways
      are with relevant genes compared to randomized version of the
      original network

Input options basic parameters
-----

     -K {integer}                
      Gene exceptions (only used for INES)

     -L<n> {integer} 
      The case (L) exceptions for the n-th matrix

     -strategy {'INES','GLONE'}  
      The strategy that will be used to extract pathways:

      INES: Will extract maximal pathways where all but at most K nodes
            are NOT active/diff. expressed in at most L cases

      GLONE: Will extract maximal pathways where the total sum of 
             NOT-active/diff. exp. cases is at most L
              
     -algo {'ACO','GREEDY','OPTIMAL'} 
       The algorithm that will be used to extract the pathways:

       GREEDY: A simple greedy algorithm, performs fast

       ACO: Ant-colony-optimization algorithm. Convergence times
               can vary depending on input size and parameters 
               (see advanced options section).

       OPTIMAL: An exact fixed-parameter-tracktability algorithm. 
               Running times increase exponentially with input 
               size(GLONE) and parameter K (INES). 
               Note: due to large running times, this will only extract 
               the best pathway

     -program {'KPM','SP','KPM_SP'}
      The program to execute: 
      KPM - KeyPathwayMiner,
      SP - Shortest paths, 
      KPM_SP - KPM followed by SP of the resulting pathways

      -sp {integer}
       Maximum number of shortest paths that will
       be reported for each pathway

      -spPathways {integer}
       Maximum number of pathways resulting from a KPM run
       that the shortest path will be executing on
       
      -spMinLength {integer}
       Minimum length of reported shortest paths

Output parameters
-----

     -fileExt {string}                
      Default file extension for output files

     -summaryFile {string}                
      Path to the file where summary will be written to

     -pathwaysFile {string}                
      Path to the file containing ALL pathways

     -pathwaysStatsFile {string}                
      Path to the file containing pathways stats

     -geneStatsFile {string}                
      Path to the file containing gene stats

     -statsFile {string}                
      Path to general stats file

     -dataStatsFile {string}                
      Path to the file where stats of the datasets will be output to

     -resultsDir {string}                
      Folder to output the result files

     -spStatsFile {string}                
      Path to shortest path stats file

     -spFile {string}                
      Path to shortest path file

     -spNodeStatsFile {string}                
      Path to shortest path node stats file

     -spEdgeStatsFile {string}                
      Path to shortest path edge stats file

     -pSingleFile {boolean}                
      If all pathways should be written to a single file or in separate ones

     -gSummary {boolean}                
      If the summary file should be generated

     -gPathways {boolean}                
      If the pathways stats file should be generated

     -gPathwayStats {boolean}                
      If the pathway stats file should be generated

     -gGeneStats {boolean}                
      If the gene stats file should be generated

     -gDataStats {boolean}                
      If the datasets stats should be generated

     -gSPStats {boolean}                
      If the shortest paths stats should be generated

     -gSPFiles {boolean}                
      If shortest paths files should be generated

     -gSPNodes {boolean}                
      If the shortest paths nodes stats should be generated

     -gSPEdges {boolean}                
      If the shortest paths edges stats should be generated

     -pGraphStats {boolean}                
      If dataset stats should be output to terminal

     -pDataStats {boolean}                
      If graph stats should be output to terminal

     -suffix {string}                
      Suffix to add to end of output files

     -runID {string}                
      Run ID

	  
Advanced options
-----

     -numProc  {integer}
      Number of threads to use (for parallel computing)

     -nodeHeuristic {'TOTAL','AVERAGE'}
      The heuristic value for each node when searching for solutions.
      This can be: 
      AVERAGE (average differentially expressed cases) or
      TOTAl (total number of differentially expressed cases)

     -combineOp {OR,AND,CUSTOM}
       How to combine multiple matrices using boolean operators.
       If CUSTOM is chosen then the logical predicate
       defined in the kpm.properties file will be used

     -combineFormula {string}
      The boolean formula used to combine the different datasets. Used
      only if combineOp == CUSTOM

      Valid operators: 
      && = AND, 
      || = OR, 
      ! = negation, 
      () = parenthesis

     -eval {boolean}
      Determines whether certain evaluation routines should run. Enabling only
      yields some statistics, Has no effect on a "normal" algorithm run other
      than slowing it down

     -maxSolutions {integer}
      Maximum number of reported pathways 
      (default is 20). Ignored if
      OPTIMAL algorithm is selected

     -doubleSolutions {boolean}
      Whether the solution array is allowed to yield multiple entries of the
      same solution 
    
     -removeBens {flag}
      If set  border exception nodes will be removed

### Advanced options ACO
     -alpha {double}
      Parameter to control the importance given 
      to the pheromone.

     -beta {double}
      Parameter to control the importance given to the
      heuristic value of the node

     -rho {double}
      Parameter that controls the pheromone decay rate

     -tauMin {double}
      Minimum pheromone that can be on a node

     -iterations {integer} 
      Maximum number of iterations

     -maxrunswithoutchange {integer}
      Maximum number of iterations allowed without
      improvement in the best solution.

     -iterationbased {boolean}
      If iteration or global best ACO strategy should be used 
      (only for GLONE). 
     
     -tradeoff {'multiplicative','additive'}
      Defines the tradeoff between pheromones and fitness when an ant picks a 
      new vertex. 
      If multiplicative, then tradeOff(a,b) = a^(alpha)*b^(beta). 
      If additive, tradeOff(a,b) = alpha * a + beta * b.

     -seed {long}
      Seed

     -solutions {integer}
      Number of solutions per iteration

     -startNodes {integer}
      Number of startnodes
      
     -rhoDecay {'CONSTANT', 'LINEAR', 'QUADRATIC', 'EXPONENTIAL'}
      The function that determines how fast the pheromone Rho should
      decay: CONSTANT, LINEAR, QUADRATIC or EXPONENTIAL

     -iterationBased {boolean}
      If TRUE uses and iteration bases ACO, if FALSE uses a global based ACO

     -localSearch {'GREEDY1', 'GREEDY2', 'OPTIMAL', 'OFF'}
      Which local search method is used to improve the results:
      GREEDY1, GREEDY2, OPTIMAL or OFF
    
     -maxRunsWithoutChange {integer}
      Maximum number of iterations allowed without improvement

### Advanced options for network robustness 
     -perturbation {quadruple of integers}
      Quadruple with parameters in this order:
        startPercent: Perturbation percentage range lower value
        stepPercent: Perturbation percentage step size
        maxPercent:  Perturbation percentage range upper value
        graphsPerStep: Number of random graphs to be created (permutations)
      Example: -perturbation=10,10,20,1
      

     -perturbation_technique {'edgeremove','edgerewire','nodeswap','noderemove'}
      Perturbation technique:
        1. edgeremove 
        2. edgerewire
        3. nodeswap
        4. noderemove

### Advanced options for batch runs
     -batch {flag}
      If set a batch run will be performed
      
     -K_batch {tripel of integers}
        Ranged gene exceptions s(only used for INES)
        MIN_K: Integer, starting value of k range or k value if k is not ranged
        INC_K: Integer, how k should be increased within the range
        MAX_K: Integer, the maximum k value, i.e. the upper limit of the range
      Example: -K_batch=1,2,3 

        
     -L<n>_batch {triple of integers}
        The ranged case (L) exceptions for the n-th matrix
        MIN_L: Integer, starting value of l range or l value if l is not ranged
        INC_L: Integer, how l should be increased within the range
        MAX_L: Integer, the maximum l value, i.e. the upper limit of the range
     Example:  -L1_batch=1,2,3 -L2_batch=1,2,3
