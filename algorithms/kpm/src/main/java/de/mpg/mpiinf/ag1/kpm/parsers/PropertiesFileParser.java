package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.Program;
import de.mpg.mpiinf.ag1.kpm.main.Main;
import de.mpg.mpiinf.ag1.kpm.utils.OutputSettings;
import de.mpg.mpiinf.ag1.kpm.utils.Separator;

import dk.sdu.kpm.Algo;
import dk.sdu.kpm.Combine;
import dk.sdu.kpm.Heuristic;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.algo.glone.LocalSearch;
import dk.sdu.kpm.algo.glone.RhoDecay;

public class PropertiesFileParser {

    private volatile KPMSettings kpmSettings;
    private String propertiesFile;

    public PropertiesFileParser(KPMSettings settings) {
        this.kpmSettings = settings;
    }

    public PropertiesFileParser(KPMSettings kpmSettings, String propertiesFile) {
        this.kpmSettings = kpmSettings;
        this.propertiesFile = propertiesFile;
    }

    public Parameters parse(String datasetFolder) {
        Parameters params = new Parameters();
        try {
            InputStream is = new FileInputStream(propertiesFile);
            //Properties props = new Properties(setDefaultProperties());
            Properties props = new Properties();
            props.load(is);
            is.close();

            /* ---- INPUT FILES ---- */
            params.GRAPH_FILE = datasetFolder + props.getProperty("graph_file");
            params.IS_GRAPH_DIRECTED = Boolean.parseBoolean(props.getProperty("is_graph_directed"));
            params.GRAPH_FILE_HAS_HEADER = Boolean.parseBoolean(props.getProperty("graph_file_has_header"));
            params.GRAPH_FILE_SEPARATOR = Separator.valueOf(props.getProperty("graph_file_separator"));

            params.DATASETS_FILE = datasetFolder + props.getProperty("datasets_file");
            params.DATASETS_FILE_HAS_HEADER = Boolean.parseBoolean(props.getProperty("datasets_file_has_header"));
            params.DATASETS_FILE_SEPARATOR = Separator.valueOf(props.getProperty("datasets_file_separator"));

            params.POSITIVE_FILE = datasetFolder + props.getProperty("positive_file");
            params.NEGATIVE_FILE = datasetFolder + props.getProperty("negative_file");

            params.MATRIX_FILES_HAVE_HEADER = Boolean.parseBoolean(props.getProperty("matrix_files_have_header"));
            params.MATRIX_FILES_SEPARATOR = Separator.valueOf(props.getProperty("matrix_files_separator"));

            params.VALIDATION_FILE = datasetFolder + props.getProperty("validation_file");

            /* ---- KPM OUTPUT FILES ---- */
            params.FILE_EXTENSION = props.getProperty("file_extension");
            params.RESULTS_FOLDER = props.getProperty("results_folder");
            params.PATHWAYS_IN_SINGLE_FILE = Boolean.parseBoolean(props.getProperty("pathways_single_file"));
            params.SUFFIX = props.getProperty("suffix");

            OutputSettings.SUMMARY_FILE = props.getProperty("summary_file");
            OutputSettings.PATHWAYS_FILE = props.getProperty("pathways_file");
            OutputSettings.PATHWAYS_STATS_FILE = props.getProperty("pathways_stats_file");
            OutputSettings.GENE_STATS_FILE = props.getProperty("gene_stats_file");
            OutputSettings.DATASETS_STATS_FILE = props.getProperty("datasets_stats_file");
            OutputSettings.GENERATE_SUMMARY_FILE = Boolean.parseBoolean(props.getProperty("generate_summary_file"));
            OutputSettings.GENERATE_PATHWAYS_FILE = Boolean.parseBoolean(props.getProperty("generate_pathways_file"));
            OutputSettings.GENERATE_PATHWAYS_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_pathways_stats_file"));
            OutputSettings.GENERATE_GENE_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_gene_stats_file"));
            OutputSettings.GENERATE_DATASETS_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_datasets_stats_file"));

            /* ---- SHORTEST PATHS OUTPUT FILES ---- */

            params.SHORTEST_PATHS_STATS_FILE = props.getProperty("shortest_paths_stats_file");

            params.SHORTEST_PATH_FILE = props.getProperty("shortest_path_file");

            params.SHORTEST_PATHS_NODE_STATS_FILE = props.getProperty("shortest_paths_node_stats_file");

            params.SHORTEST_PATHS_EDGE_STATS_FILE = props.getProperty("shortest_paths_edge_stats_file");

            params.GENERATE_SHORTEST_PATHS_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_shortest_paths_stats_file"));

            params.GENERATE_SHORTEST_PATH_FILES = Boolean.parseBoolean(props.getProperty("generate_shortest_paths_files"));

            params.GENERATE_SHORTEST_PATHS_NODE_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_shortest_paths_node_stats_file"));

            params.GENERATE_SHORTEST_PATHS_EDGE_STATS_FILE = Boolean.parseBoolean(props.getProperty("generate_shortest_paths_edge_stats_file"));

            /* ---- BASIC params -------*/
            kpmSettings.GENE_EXCEPTIONS = Integer.parseInt(props.getProperty("gene_exceptions"));
            String strategy = props.getProperty("strategy");
            String algorithm = props.getProperty("algorithm");
            if (strategy.equals("INES")) {
                if (algorithm.equals("GREEDY")) {
                    kpmSettings.ALGO = Algo.GREEDY;
                } else if (algorithm.equals("ACO")) {
                    kpmSettings.ALGO = Algo.LCG;
                } else if (algorithm.equals("OPTIMAL")) {
                    kpmSettings.ALGO = Algo.OPTIMAL;
                } else {
                    kpmSettings.ALGO = Algo.GREEDY;
                }
            } else if (strategy.equals("GLONE")) {
                if (algorithm.equals("GREEDY")) {
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMGREEDY;
                } else if (algorithm.equals("ACO")) {
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMACO;
                } else if (algorithm.equals("OPTIMAL")) {
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMOPTIMAL;
                } else {
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMGREEDY;
                }

            } else {
                if (algorithm.equals("GREEDY")) {
                    kpmSettings.ALGO = Algo.GREEDY;
                } else if (algorithm.equals("ACO")) {
                    kpmSettings.ALGO = Algo.LCG;
                } else if (algorithm.equals("OPTIMAL")) {
                    kpmSettings.ALGO = Algo.OPTIMAL;
                } else {
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMGREEDY;
                }
            }

            params.PROGRAM = Program.valueOf(props.getProperty("program"));

            /* ------ SHORTEST PATHS params ---------------- */
            params.PATHWAYS_SHORTEST_PATHS =
                    Integer.parseInt(props.getProperty("pathways_shortest_paths"));

            params.SHORTEST_PATHS_REPORTED =
                    Integer.parseInt(props.getProperty("shortest_paths_reported"));

            params.MINIMUM_LENGTH_SHORTEST_PATHS =
                    Integer.parseInt(props.getProperty("minimum_length_shortest_paths"));

            params.SORT_SHORTEST_PATHS_BY =
                    de.mpg.mpiinf.ag1.kpm.shortestpath.SortShortestPaths.valueOf(props.getProperty("sort_shortest_paths_by"));
            /* ------ ADVANCED params (GENERAL) ------------ */

            int max = Runtime.getRuntime().availableProcessors();
            String procProp = props.getProperty("processors");
            if (procProp.equals("MAX")) {
                kpmSettings.NUMBER_OF_PROCESSORS = max;
            } else {
                int val = Integer.parseInt(procProp);
                if (val > max) {
                    kpmSettings.NUMBER_OF_PROCESSORS = max;
                } else if (val < 1) {
                    kpmSettings.NUMBER_OF_PROCESSORS = 1;
                } else {
                    kpmSettings.NUMBER_OF_PROCESSORS = val;
                }
            }
            kpmSettings.NODE_HEURISTIC_VALUE =
                    Heuristic.valueOf(props.getProperty("node_heuristic"));
            kpmSettings.EVAL =
                    Boolean.parseBoolean(props.getProperty("eval"));
            kpmSettings.DOUBLE_SOLUTIONS_ALLOWED =
                    Boolean.parseBoolean(props.getProperty("double_solutions"));
            kpmSettings.COMBINE_OPERATOR =
                    Combine.valueOf(props.getProperty("combine_operator"));
            kpmSettings.COMBINE_FORMULA =
                    props.getProperty("combine_formula");

            /* ------ ADVANCED PARAMETERS (ONLY FOR ACO ALGORITHM) ------------ */
            kpmSettings.ALPHA = Double.parseDouble(props.getProperty("alpha"));
            kpmSettings.BETA = Double.parseDouble(props.getProperty("beta"));
            kpmSettings.RHO = Double.parseDouble(props.getProperty("rho"));
            kpmSettings.RHO_DECAY = RhoDecay.valueOf(props.getProperty("rho_decay"));
            kpmSettings.TAU_MIN = Double.parseDouble(props.getProperty("tau_min"));
            kpmSettings.MULTIPLICATIVE_TRADEOFF =
                    Boolean.parseBoolean(props.getProperty("multiplicative_tradeoff"));
            kpmSettings.ITERATION_BASED =
                    Boolean.parseBoolean(props.getProperty("iteration_based"));
            int maxIterations = Integer.parseInt(props.getProperty("max_iterations"));
            if (maxIterations < 0) {
                maxIterations = Integer.MAX_VALUE;
            }
            kpmSettings.MAX_ITERATIONS = maxIterations;
            int maxIterWOChange = Integer.parseInt(props.getProperty("max_iterations_wo_change"));
            if (maxIterWOChange < 0) {
                maxIterWOChange = Integer.MAX_VALUE;
            }
            kpmSettings.MAX_RUNS_WITHOUT_CHANGE = maxIterWOChange;
            kpmSettings.NUMBER_OF_SOLUTIONS_PER_ITERATION = Integer.parseInt(props.getProperty("num_solutions_per_iteration"));
            kpmSettings.NUM_STARTNODES = Integer.parseInt(props.getProperty("num_startnodes"));
            kpmSettings.L_SEARCH = LocalSearch.valueOf(props.getProperty("local_search"));


        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException io) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, io);
        }

        return params;
    }
}
