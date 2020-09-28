package de.mpg.mpiinf.ag1.kpm;

import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.KPMGraph;

public class KPMParameterConfigurator {
    private volatile KPMGraph graph;

    private volatile KPMAlgorithm algorithm;

    private volatile KPMStrategy strategy;

    private String prefix = "Configuration is incorrect:\t";

    private volatile KPMSettings kpmSettings;

    public KPMParameterConfigurator(Parameters params, KPMSettings settings) {
        this.kpmSettings = settings;
        this.graph = kpmSettings.MAIN_GRAPH;
        this.algorithm = params.ALGORITHM;
        this.strategy = params.STRATEGY;
    }

    public boolean areSettingsValid() {
        kpmSettings.STARTING_TIME = System.nanoTime();

        if (graph == null) {
            System.out.println(prefix + "No graph loaded.");
            return false;
        }

        int maxNumberOfGeneExceptions = graph.getNodeIdSet().size();

        if (kpmSettings.IS_BATCH_RUN) {
            if (strategy == KPMStrategy.INES) {
                if (kpmSettings.MAX_K <= 0) {
                    System.out.println(prefix + "Invalid MAX_K.");
                    return false;
                }

                if (kpmSettings.MIN_K > kpmSettings.MAX_K) {
                    System.out.println(prefix + "Invalid MIN_K. It is larger than MAX_K.");
                    return false;
                }

                if (kpmSettings.INC_K == 0) {
                    System.out.println(prefix + "Invalid INC_K. Must be larger than 0.");
                    return false;
                }

                if (kpmSettings.MAX_K - kpmSettings.INC_K < kpmSettings.MIN_K) {
                    System.out.println(prefix + "Invalid INC_K. Incrementation must be in range.");
                    return false;
                }
            } else {
                kpmSettings.MIN_K = kpmSettings.MAX_K = kpmSettings.INC_K = 0;
                kpmSettings.GENE_EXCEPTIONS = 0;
            }
            // Checking L's
            for (String internalId : kpmSettings.MIN_L.keySet()) {
                if (kpmSettings.MAX_L.get(internalId) <= 0) {
                    System.out.println(prefix + "Invalid MAX_L.");
                    return false;
                }
                if (kpmSettings.MIN_L.get(internalId) > kpmSettings.MAX_L.get(internalId)) {
                    System.out.println(prefix + "Invalid MIN_L. It is larger than MAX_L.");
                    return false;
                }

                if (kpmSettings.INC_L.get(internalId) == 0) {
                    System.out.println(prefix + "Invalid INC_L. Must be larger than 0.");
                    return false;
                }

                if ((kpmSettings.MAX_L.get(internalId) - kpmSettings.INC_L.get(internalId) < kpmSettings.MIN_L.get(internalId))
                        && !(kpmSettings.MAX_L.get(internalId).equals(kpmSettings.INC_L.get(internalId))
                        && kpmSettings.MAX_L.get(internalId).equals(kpmSettings.MIN_L.get(internalId)))) {
                    // IN case INC_L is not in range and L parameters are not equals
                    System.out.println(prefix + "Invalid INC_L. Incrementation must be in range.");
                    return false;
                }
            }
        } else {
            if (strategy == KPMStrategy.INES) {
                if (maxNumberOfGeneExceptions < kpmSettings.GENE_EXCEPTIONS) {
                    System.out.println(prefix + "Invalid number of GENE_EXCEPTIONS. Must be less than the number of nodes.");
                    return false;
                }
            }
        }

        return true;
    }
}
