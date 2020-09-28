package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
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
import dk.sdu.kpm.perturbation.PerturbationService;
import dk.sdu.kpm.perturbation.IPerturbation.PerturbationTags;

public class ArgsParametersParser {

    private volatile KPMSettings kpmSettings;

    //List of all algorithms to test if the right algorithm was provided
    private final List<String> algoList = Arrays.asList("GREEDY", "ACO", "OPTIMAL");
    private String algorithm = "GREEDY";

    //List of all strategies to test if the right strategy was provided
    private final List<String> strategyList = Arrays.asList("INES", "GLONE");
    private String strategy = "INES";

    private HashMap<String, String> id2path = new HashMap<String, String>();
    private HashMap<String, Integer> id2param = new HashMap<String, Integer>();


    public ArgsParametersParser(KPMSettings settings) {
        this.kpmSettings = settings;
    }

    /**
     * Parses command line arguments and
     * saves them to the corresponding objects.
     *
     * @param args   Arguments passed from the commandline.
     * @param params Object with parameters for the standalone execution.
     * @return Parameters object.
     * @throws Exception
     */
    public Parameters parse(String[] args, Parameters params) throws Exception {
        //Read command line arguments
        parseParameters(args, params);

        if (params.PROGRAM == Program.SP) {
            // do nothing
        } else if (!id2path.isEmpty()) {
            for (String id : id2path.keySet()) {
                String file = id2path.get(id);
                System.out.println(id + ": " + file);

                if (!(new File(file)).isFile()) {
                    System.out.println(file + " not found !");
                    System.exit(-1);
                }

            }
            kpmSettings.MATRIX_FILES_MAP = id2path;

            if (kpmSettings.IS_BATCH_RUN) {
                List<String> invalids = new ArrayList<String>();
                for (String id : id2param.keySet()) {
                    if (id.startsWith("L")) {
                        invalids.add(id);
                    }
                }
                for (String invalidID : invalids) {
                    id2param.remove(invalidID);
                }
            }
            kpmSettings.CASE_EXCEPTIONS_MAP = id2param;
        } else if (params.DATASETS_FILE == null) {
            System.out.println("No datasets file was specified.");
            System.exit(-1);
        } else if (!new File(params.DATASETS_FILE).isFile()) {
            System.out.println("Datasets file " + params.DATASETS_FILE + " does not exist.");
            System.exit(-1);
        } else {
            DatasetsFileParser dfp = new DatasetsFileParser(kpmSettings);
            params = dfp.parse(params.DATASETS_FILE_SEPARATOR.charValue(), params.DATASETS_FILE_HAS_HEADER, params);
        }

        //Compares the number of matrices against the number of given L parameters
        if (kpmSettings.MIN_L.keySet().size() != kpmSettings.MATRIX_FILES_MAP.size()) {
            System.out.println(String.format(
                    "\nThe were found setup for %d L-parameters, this amount does not match the %d found for the" +
                            " matrix files map.\n\nKPM will now terminate.",
                    kpmSettings.MIN_L.keySet().size(),
                    kpmSettings.MATRIX_FILES_MAP.size()
            ));
            System.exit(-1);
        }

        //Set the right algorithm depending on the strategy
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

        // Custom combine formula for combining the different datasets
        if (kpmSettings.COMBINE_OPERATOR.toString().equals("CUSTOM")) {
            System.out.println("Combine Formula : " + kpmSettings.COMBINE_FORMULA);
        }

        return params;
    }

    private void parseParameters(String[] args, Parameters params) throws Exception {
        for (String arg : args) {
            // Split argument on equals symbol
            String[] options = arg.split("=");

            // Input files
            if (options[0].equals("-graphFile")) {
                params.GRAPH_FILE = options[1];
            } else if (options[0].startsWith("-isDirected")) {
                params.IS_GRAPH_DIRECTED = Boolean.parseBoolean(options[1]);
            } else if (options[0].startsWith("-matrix")) {
                String id = "L" + options[0].substring(7);
                String path = options[1];
                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);
                id2path.put(internalID, path);
            } else if (options[0].equals("-gfHeader")) {
                params.GRAPH_FILE_HAS_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gfSep")) {
                params.GRAPH_FILE_SEPARATOR = Separator.valueOf(options[1]);
            } else if (options[0].equals("-datasetsFile")) {
                params.DATASETS_FILE = options[1];
            } else if (options[0].equals("-dfHeader")) {
                params.DATASETS_FILE_HAS_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-dfSep")) {
                params.DATASETS_FILE_SEPARATOR = Separator.valueOf(options[1]);
            } else if (options[0].equals("-positiveFile")) {
                params.POSITIVE_FILE = options[1];
            } else if (options[0].equals("-negativeFile")) {
                params.NEGATIVE_FILE = options[1];
            } else if (options[0].equals("-mfHeader")) {
                params.MATRIX_FILES_HAVE_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-mfSep")) {
                params.MATRIX_FILES_SEPARATOR = Separator.valueOf(options[1]);
                // OUTPUT FILES
            } else if (options[0].equals("-fileExt")) {
                params.FILE_EXTENSION = options[1];
            } else if (options[0].equals("-summaryFile")) {
                OutputSettings.SUMMARY_FILE = options[1];
            } else if (options[0].equals("-pathwaysFile")) {
                OutputSettings.PATHWAYS_FILE = options[1];
            } else if (options[0].equals("-pathwaysStatsFile")) {
                OutputSettings.PATHWAYS_STATS_FILE = options[1];
            } else if (options[0].equals("-geneStatsFile")) {
                OutputSettings.GENE_STATS_FILE = options[1];
            } else if (options[0].equals("statsFile")) {
                OutputSettings.GENERAL_STATS_FILE = options[1];
            } else if (options[0].equals("-dataStatsFile")) {
                OutputSettings.DATASETS_STATS_FILE = options[1];
            } else if (options[0].equals("-resultsDir")) {
                params.RESULTS_FOLDER = options[1];
            } else if (options[0].equals("-spStatsFile")) {
                params.SHORTEST_PATHS_STATS_FILE = options[1];
            } else if (options[0].equals("-spFile")) {
                params.SHORTEST_PATH_FILE = options[1];
            } else if (options[0].equals("-spNodeStatsFile")) {
                params.SHORTEST_PATHS_NODE_STATS_FILE = options[1];
            } else if (options[0].equals("-spEdgeStatsFile")) {
                params.SHORTEST_PATHS_EDGE_STATS_FILE = options[1];
            } else if (options[0].equals("-pSingleFile")) {
                params.PATHWAYS_IN_SINGLE_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSummary")) {
                OutputSettings.GENERATE_SUMMARY_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gPathways")) {
                OutputSettings.GENERATE_PATHWAYS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gPathwayStats")) {
                OutputSettings.GENERATE_PATHWAYS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gGeneStats")) {
                OutputSettings.GENERATE_GENE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gDataStats")) {
                OutputSettings.GENERATE_DATASETS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPStats")) {
                params.GENERATE_SHORTEST_PATHS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPFiles")) {
                params.GENERATE_SHORTEST_PATH_FILES = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPNodes")) {
                params.GENERATE_SHORTEST_PATHS_NODE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPEdges")) {
                params.GENERATE_SHORTEST_PATHS_EDGE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("runID")) {
                params.RUN_ID = options[1];
                // OUTPUT TO TERMINAL
            } else if (options[0].equals("-pGraphStats")) {
                params.PRINT_GRAPH_STATS = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-pDataStats")) {
                params.PRINT_DATASETS_STATS = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-suffix")) {
                params.SUFFIX = options[1];
                // BASIC params
            } else if (options[0].equals("-program")) {
                params.PROGRAM = Program.valueOf(options[1]);
            } else if (options[0].equals("-sp")) {
                params.SHORTEST_PATHS_REPORTED = Integer.valueOf(options[1]);
            } else if (options[0].equals("-spPathways")) {
                params.PATHWAYS_SHORTEST_PATHS = Integer.valueOf(options[1]);
            } else if (options[0].equals("-spMinLength")) {
                params.MINIMUM_LENGTH_SHORTEST_PATHS = Integer.valueOf(options[1]);
            } else if (options[0].equals("-K")) {
                kpmSettings.GENE_EXCEPTIONS = Integer.parseInt(options[1]);
            } else if (options[0].startsWith("-L") && !options[0].contains("batch")) {
                String id = "L" + options[0].substring(2);
                int l = Integer.parseInt(options[1]);
                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);
                kpmSettings.MIN_L.put(internalID, l);
                kpmSettings.INC_L.put(internalID, l);
                kpmSettings.MAX_L.put(internalID, l);
                id2param.put(internalID, l);
            } else if (options[0].equals("-algo")) {
                if (!algoList.contains(options[1])) {
                    System.err.println(options[1]
                            + " is not a valid algorithm. " + ("Valid options are: GREEDY, ACO or OPTIMAL"));
                    System.exit(-1);
                } else {
                    algorithm = options[1];
                }

            } else if (options[0].equals("-strategy")) {
                if (!strategyList.contains(options[1])) {
                    System.err.println(options[1]
                            + " is not a valid strategy. " + ("Valid options are: INES or GLONE"));
                    System.exit(-1);

                } else {
                    strategy = options[1];
                }

                // ADVANCED PARAMETERS (GENERAL)

            } else if (options[0].equals("-numProc")) {
                int input = Integer.parseInt(options[1]);
                int available = Runtime.getRuntime().availableProcessors();
                if (input < 1) {
                    kpmSettings.NUMBER_OF_PROCESSORS = 1;
                } else if (input > available) {
                    kpmSettings.NUMBER_OF_PROCESSORS = available;
                } else {
                    kpmSettings.NUMBER_OF_PROCESSORS = input;
                }
            } else if (options[0].equals("-nodeHeuristic")) {
                kpmSettings.NODE_HEURISTIC_VALUE = Heuristic.valueOf(options[1]);
            } else if (options[0].equals("-combineOp")) {
                kpmSettings.COMBINE_OPERATOR = Combine.valueOf(options[1]);
            } else if (options[0].equals("-combineFormula")) {
                kpmSettings.COMBINE_FORMULA = options[1];
            } else if (options[0].equals("-eval")) {
                kpmSettings.EVAL = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-maxSolutions")) {
                kpmSettings.NUM_SOLUTIONS = Integer.parseInt(options[1]);
            } else if (options[0].equals("-doubleSolutions")) {
                kpmSettings.DOUBLE_SOLUTIONS_ALLOWED = Boolean.parseBoolean(options[1]);

                // ADVANCED PARAMETERS (ACO ONLY)
            } else if (options[0].equals("-alpha")) {
                kpmSettings.ALPHA = Double.parseDouble(options[1]);
            } else if (options[0].equals("-beta")) {
                kpmSettings.BETA = Double.parseDouble(options[1]);
            } else if (options[0].equals("-rho")) {
                kpmSettings.RHO = Double.parseDouble(options[1]);
            } else if (options[0].equals("-iterations")) {
                if (Integer.parseInt(options[1]) == 0) {
                    kpmSettings.MAX_ITERATIONS = Integer.MAX_VALUE;
                } else {
                    kpmSettings.MAX_ITERATIONS = Integer.parseInt(options[1]);
                }
            } else if (options[0].equals("-seed")) {
                kpmSettings.SEED = Long.parseLong(options[1]);
                kpmSettings.R = new Random(kpmSettings.SEED);
            } else if (options[0].equals("-tradeoff")
                    && options[1].equals("multiplicative")) {
                kpmSettings.MULTIPLICATIVE_TRADEOFF = true;
            } else if (options[0].equals("-tradeoff")
                    && options[1].equals("additive")) {
                kpmSettings.MULTIPLICATIVE_TRADEOFF = false;
            } else if (options[0].equals("-solutions")) {
                kpmSettings.NUMBER_OF_SOLUTIONS_PER_ITERATION = Integer.parseInt(options[1]);
            } else if (options[0].equals("-startNodes")) {
                kpmSettings.NUM_STARTNODES = Integer.parseInt(options[1]);
            } else if (options[0].equals("-rhoDecay")) {
                try {
                    kpmSettings.RHO_DECAY = RhoDecay.valueOf(options[1]);
                } catch (IllegalArgumentException e) {
                    System.err.println(options[1]
                            + " is not a valid rho decay function. "
                            + Algo.values());
                    System.exit(-1);
                }
            } else if (options[0].equals("-iterationBased")) {
                try {
                    kpmSettings.ITERATION_BASED = Boolean.parseBoolean(options[1]);
                } catch (IllegalArgumentException e) {
                    System.err.println(options[1]
                            + " is not a valid local search method. "
                            + LocalSearch.values());
                    System.exit(-1);
                }
            } else if (options[0].equals("-localSearch")) {
                kpmSettings.L_SEARCH = LocalSearch.valueOf(options[1]);
            } else if (options[0].equals("-maxRunsWithoutChange")) {
                if (Integer.parseInt(options[1]) == 0) {
                    kpmSettings.MAX_RUNS_WITHOUT_CHANGE = Integer.MAX_VALUE;
                } else {
                    kpmSettings.MAX_RUNS_WITHOUT_CHANGE = Integer.parseInt(options[1]);
                }
            } else if (options[0].equals("-tauMin")) {
                kpmSettings.TAU_MIN = Double.parseDouble(options[1]);
            } else if (options[0].equals("-batch")) {
                System.out.println("Should be batch run");
                kpmSettings.IS_BATCH_RUN = true;
            } else if (options[0].equals("-removeBens")) {
                kpmSettings.REMOVE_BENs = true;
            } else if (options[0].equals("-validationFile")) {
                params.VALIDATION_FILE = options[1];
            } else if (options[0].equals("-perturbation")) {
                params.IS_PERTURBATION_RUN = true;
                String[] values = options[1].split(",");

                if (values.length != 4) {
                    throw new Exception("Invalid settings for " + options[0]);
                }
                params.MIN_PERCENTAGE = Integer.parseInt(values[0]);
                params.STEP_PERCENTAGE = Integer.parseInt(values[1]);
                params.MAX_PERCENTAGE = Integer.parseInt(values[2]);
                params.GRAPHS_PER_STEP = Integer.parseInt(values[3]);
            }
            // Added setting for batch K
            else if (options[0].equals("-K_batch")) {
                String[] values = options[1].split(",");

                if (values.length != 3) {
                    throw new Exception("Invalid settings for " + options[0]);
                }

                kpmSettings.MIN_K = Integer.parseInt(values[0]);
                kpmSettings.INC_K = Integer.parseInt(values[1]);
                kpmSettings.MAX_K = Integer.parseInt(values[2]);

            } else if (options[0].matches("-L[1-9][0-9]*[_]batch")) {
                // If batch option is set assign ranged L value to n th matrix
                String id = options[0].substring(1, options[0].indexOf('_')); // get number of nth L parmeter
                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);
                String[] values = options[1].split(",");

                if (values.length != 3) {
                    throw new Exception("Invalid settings for " + options[0]);
                }

                kpmSettings.MIN_L.put(internalID, Integer.parseInt(values[0]));
                kpmSettings.INC_L.put(internalID, Integer.parseInt(values[1]));
                kpmSettings.MAX_L.put(internalID, Integer.parseInt(values[2]));
                kpmSettings.VARYING_L_ID.add(internalID);


            } else if (options[0].equals("-perturbation_technique")) {
                if (options[1].equals("edgeremove")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.EdgeRemoval);

                } else if (options[1].equals("edgerewire")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.EdgeRewire);

                } else if (options[1].equals("nodeswap")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.NodeSwap);

                } else if (options[1].equals("noderemove")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.NodeRemoval);

                }

                if (params.PERTURBATION == null) {
                    throw new Exception("Invalid perturbation technique.");
                }

                System.out.println("Perturbation technique: " + params.PERTURBATION.getName());
            } else if (options[0].equals("-help")) {
                printHelp();
                System.exit(0);
            } else {
                System.out.println("Unknown command line argument: "
                        + options[0]);
                System.out.println("Please type \"-help\" for a list of commands and usage.");
                System.exit(-1);
            }
        }
    }

    private static void printHelp() {
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader("README.md"));
            String line = "";
            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ioe);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
