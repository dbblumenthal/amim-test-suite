package de.mpg.mpiinf.ag1.kpm;

import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.utils.OutputSettings;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.Result;

public class ObsoleteMethods {
    private static Properties setDefaultProperties(Parameters params, KPMSettings kpmSettings) {
        Properties def = new Properties();
        /* ----- INPUT FILES DEFAULTS---------*/
        def.put("graph_file_has_header", String.valueOf(params.GRAPH_FILE_HAS_HEADER));
        def.put("graph_file_separator", params.GRAPH_FILE_SEPARATOR.toString());
        def.put("datasets_file_has_header", String.valueOf(params.DATASETS_FILE_HAS_HEADER));
        def.put("datasets_file_separator", params.DATASETS_FILE_SEPARATOR.toString());
        def.put("matrix_files_have_header", String.valueOf(params.MATRIX_FILES_HAVE_HEADER));
        def.put("matrix_files_separator", params.MATRIX_FILES_SEPARATOR.toString());

        /* ----- OUTPUT FILES DEFAULTS---------*/
        OutputSettings.SUMMARY_FILE = params.RESULTS_FOLDER + params.FILE_SEPARATOR +OutputSettings.SUMMARY_FILE;
        def.put("summary_file",OutputSettings.SUMMARY_FILE);
        OutputSettings.PATHWAYS_FILE = params.RESULTS_FOLDER + params.FILE_SEPARATOR +OutputSettings.PATHWAYS_FILE;
        def.put("pathways_file",OutputSettings.PATHWAYS_FILE);
        OutputSettings.PATHWAYS_STATS_FILE = params.RESULTS_FOLDER + params.FILE_SEPARATOR +OutputSettings.PATHWAYS_STATS_FILE;
        def.put("pathways_stats_file",OutputSettings.PATHWAYS_STATS_FILE);
        OutputSettings.GENE_STATS_FILE = params.RESULTS_FOLDER + params.FILE_SEPARATOR +OutputSettings.GENE_STATS_FILE;
        def.put("gene_stats_file",OutputSettings.GENE_STATS_FILE);
        OutputSettings.DATASETS_STATS_FILE = params.RESULTS_FOLDER + params.FILE_SEPARATOR + OutputSettings.DATASETS_STATS_FILE;
        def.put("datasets_stats_file", OutputSettings.DATASETS_STATS_FILE);
        def.put("pathways_single_file", params.PATHWAYS_IN_SINGLE_FILE);
        def.put("generate_summary_file", String.valueOf(OutputSettings.GENERATE_SUMMARY_FILE));
        def.put("generate_pathways_file", String.valueOf(OutputSettings.GENERATE_PATHWAYS_FILE));
        def.put("generate_pathways_stats_file", String.valueOf(OutputSettings.GENERATE_PATHWAYS_STATS_FILE));
        def.put("generate_gene_stats_file", String.valueOf(OutputSettings.GENERATE_GENE_STATS_FILE));
        def.put("generate_datasets_stats_file", String.valueOf(OutputSettings.GENERATE_DATASETS_STATS_FILE));
        def.put("suffix", params.SUFFIX);
        def.put("file_extension", params.FILE_EXTENSION);

        /* ----- BASIC params DEFAULTS -----*/
        def.put("gene_exceptions", String.valueOf(kpmSettings.GENE_EXCEPTIONS));
        def.put("strategy", "INES");
        def.put("algorithm", "GREEDY");
        def.put("max_solutions", String.valueOf(kpmSettings.NUM_SOLUTIONS));

        /* ----- ADANCED params (GENERAL) DEFAULTS ----*/
        def.put("processors", String.valueOf(kpmSettings.NUMBER_OF_PROCESSORS));
        def.put("node_heuristic", kpmSettings.NODE_HEURISTIC_VALUE.toString());
        def.put("eval", String.valueOf(kpmSettings.EVAL));
        def.put("double_solutions", String.valueOf(kpmSettings.DOUBLE_SOLUTIONS_ALLOWED));
        def.put("combine_operator", String.valueOf(kpmSettings.COMBINE_OPERATOR));
        //           def.put("combine_formula", KPMParameters.COMBINE_FORMULA);

        /* ----- ADANCED params (ACO) DEFAULTS ----*/
        def.put("alpha", String.valueOf(kpmSettings.ALPHA));
        def.put("beta", String.valueOf(kpmSettings.BETA));
        def.put("rho", String.valueOf(kpmSettings.RHO));
        def.put("rho_decay", kpmSettings.RHO_DECAY.toString());
        def.put("tau_min", String.valueOf(kpmSettings.TAU_MIN));
        def.put("multiplicative_tradeoff", String.valueOf(kpmSettings.MULTIPLICATIVE_TRADEOFF));
        def.put("iteration_based", String.valueOf(kpmSettings.ITERATION_BASED));
        def.put("max_iterations", String.valueOf(kpmSettings.MAX_ITERATIONS));
        def.put("max_iterations_wo_change", String.valueOf(kpmSettings.MAX_RUNS_WITHOUT_CHANGE));
        def.put("num_solutions_per_iterations", String.valueOf(kpmSettings.NUMBER_OF_SOLUTIONS_PER_ITERATION));
        def.put("num_startnodes", String.valueOf(kpmSettings.NUM_STARTNODES));
        def.put("local_search", kpmSettings.L_SEARCH.toString());

        return def;
    }


    /**
     * Removes double results, if there are any. More specifically, for each
     * element e1 in the list, it removes all other elements e2 where
     * e1.equals(e2) yields true.
     *
     * @param results the list where the double entries are to be removed. List
     * has to be sorted.
     * @return the list without double entries. the input list is not changed.
     */
    private static List<Result> removeDoubleSolutions(List<Result> results) {
        List<Result> toReturn = new ArrayList<Result>();
        int i = 0;

        while (i < results.size()) {
            Result r = results.get(i);
            if (!toReturn.contains(r)) {
                toReturn.add(r);
            }
            i++;
        }

        return toReturn;
    }
}