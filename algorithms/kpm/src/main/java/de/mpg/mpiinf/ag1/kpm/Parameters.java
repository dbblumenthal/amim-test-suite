package de.mpg.mpiinf.ag1.kpm;

import java.util.HashMap;
import java.util.Map;

import de.mpg.mpiinf.ag1.kpm.shortestpath.SortShortestPaths;
import de.mpg.mpiinf.ag1.kpm.utils.Parser;
import de.mpg.mpiinf.ag1.kpm.utils.Separator;
import dk.sdu.kpm.graph.KPMGraph;
import dk.sdu.kpm.perturbation.IPerturbation;

public class Parameters {
  
	public Parameters() {
	}

	/* ------ INPUT FILES ------------ */

    // Map storing the paths to the multipe expression files
    // NOTE1: The key ID's must be the same as the key ID's in
    // CASE_EXCEPTIONS_MAP
    // NOTE2: Will only be used if MULTI == True
        public volatile Map<String, String> MATRIX_FILES_MAP = new HashMap<String, String>();
	
        // Where the graph file is located
        public String GRAPH_FILE = "sampleNetwork.sif";
        
        // If graph should be considered directed
        public boolean IS_GRAPH_DIRECTED = false;
        
        // If graph file has header
        public boolean GRAPH_FILE_HAS_HEADER = false;
        
        // Separating character for columns in graph file
        public Separator GRAPH_FILE_SEPARATOR = Separator.TAB;
        
        // Where the datasets file is located
        public String DATASETS_FILE = "datasets_file.txt";
        
        // If graph file has header
        public boolean DATASETS_FILE_HAS_HEADER = true;
        
        // Separating character for columns in graph file
        public Separator DATASETS_FILE_SEPARATOR = Separator.TAB;
        
        // The file containing the genes for the positive list
        public String POSITIVE_FILE = "positive_list.txt";
        
        // The file containing the genes for the negative list
        public String NEGATIVE_FILE = "negative_list.txt";
        
        public boolean MATRIX_FILES_HAVE_HEADER = false;
        
        public Separator MATRIX_FILES_SEPARATOR = Separator.TAB;
        
        public String FILE_SEPARATOR = System.getProperty("file.separator");
        
        public String VALIDATION_FILE = "COAD-VAL-ENTREZ.txt";
        
        /* ------ KPM OUTPUT FILES ------------ */
        
        public String FILE_EXTENSION = ".txt";
        
        public String RESULTS_FOLDER = "results";
              
        public boolean PATHWAYS_IN_SINGLE_FILE = false;
       

        public boolean IS_PERTURBATION_RUN = false;
        
        public int GRAPHS_PER_STEP = 0;
        
        public int MIN_PERCENTAGE = 0;
        
        public int STEP_PERCENTAGE = 0;
        
        public int MAX_PERCENTAGE = 0;
                
        
        public String SUFFIX = "KPM";
        
        public boolean PRINT_GRAPH_STATS = true;
        
        public boolean PRINT_DATASETS_STATS = true;
                
        public String RUN_ID = "";
 
        
        /* ------ BASIC PARAMETERS ------------ */

        // The program to execute: KPM - KeyPathwayMiner, SP - Shortest paths, 
        // KPM_SP - KPM followed by SP of the resulting pathways
        public Program PROGRAM = Program.KPM;   
       
        // Default strategy
        public KPMStrategy STRATEGY = KPMStrategy.INES;
        
        // Default algorithm        
        public KPMAlgorithm ALGORITHM = KPMAlgorithm.GREEDY;
        
        public KPMGraph INPUT_GRAPH = null;

        public Parser PARSER = null;
        
        public IPerturbation PERTURBATION = null;

        /* --------- SHORTEST PATHS PARAMETERS ---------------- */
        // Maximum number of pathways resulting from a KPM run
        // that the shortest path will be executing on
        public int PATHWAYS_SHORTEST_PATHS = 1;
        
        // Maximum number of shortest paths that will
        // be reported for each pathway
        public int SHORTEST_PATHS_REPORTED = 10;
        
        // Minimum length of reported shortest paths
        public int MINIMUM_LENGTH_SHORTEST_PATHS = 2;
        
        // Sort criteria for shortest paths (descending)
        public SortShortestPaths SORT_SHORTEST_PATHS_BY = SortShortestPaths.LENGTH;
        
        /* ------ Shortest paths output files ------ */
        public String SHORTEST_PATHS_STATS_FILE = "shortest_path_stats";
        
        public String SHORTEST_PATH_FILE = "shortest_path";
        
        public String SHORTEST_PATHS_NODE_STATS_FILE =  "shortest_path_node_stats";
        
        public String SHORTEST_PATHS_EDGE_STATS_FILE =  "shortest_path_edge_stats";
        
        public boolean GENERATE_SHORTEST_PATHS_STATS_FILE = true;
        
        public boolean GENERATE_SHORTEST_PATH_FILES = true;
        
        public boolean GENERATE_SHORTEST_PATHS_NODE_STATS_FILE = true;
        
        public boolean GENERATE_SHORTEST_PATHS_EDGE_STATS_FILE =  true;
        
}


