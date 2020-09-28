package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.File;

import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.utils.OutputSettings;
import dk.sdu.kpm.Algo;
import dk.sdu.kpm.KPMSettings;

public class InputFileParser {

	private volatile KPMSettings kpmSettings;

	public InputFileParser(KPMSettings settings){
		this.kpmSettings = settings;
	}

	public Parameters parse(Parameters params) {
		boolean hasErrors = false;
		if (params.GRAPH_FILE == null) {
			System.out.println("ERROR: Please specify a graph file !");
			hasErrors = true;
		} else  {
			File graph_file = new File(params.GRAPH_FILE);
			if (!graph_file.isFile()){
				System.out.println("ERROR: The specified graph file does not exist.\nPath: "+graph_file.getAbsolutePath());
				hasErrors = true;
			}
		}

		if (hasErrors) {
			System.exit(-1);
		}

		if (kpmSettings.ALGO == Algo.EXCEPTIONSUMACO
				|| kpmSettings.ALGO == Algo.EXCEPTIONSUMGREEDY
				|| kpmSettings.ALGO == Algo.EXCEPTIONSUMOPTIMAL) {
			kpmSettings.GENE_EXCEPTIONS = 0;
		}

		if (params.POSITIVE_FILE != null) {
			if (!(new File(params.POSITIVE_FILE)).exists()) {
				System.out.println("WARNING: The specified positive gene list file does not exist.");
				params.POSITIVE_FILE = null;
			}
		}

		if (params.NEGATIVE_FILE != null) {
			if (!(new File(params.NEGATIVE_FILE)).exists()) {
				System.out.println("WARNING: The specified negative gene list file does not exist.");
				params.NEGATIVE_FILE = null;
			}
		}

		if (OutputSettings.SUMMARY_FILE != null) {
			int count = 0;
			while ((new File(OutputSettings.SUMMARY_FILE + "-"
					+ params.SUFFIX + params.FILE_EXTENSION)).exists()) {
				count++;
				String name = OutputSettings.SUMMARY_FILE;
				if (count == 1) {
					OutputSettings.SUMMARY_FILE += "(" + String.valueOf(count) + ")";
				} else {
					int aux = count - 1;
					OutputSettings.SUMMARY_FILE = name.replace("(" + aux + ")", "(" + count + ")");
				}
			}
		}

		if (OutputSettings.PATHWAYS_FILE != null) {
			int count = 0;
			while ((new File(OutputSettings.PATHWAYS_FILE + "-"
					+ params.SUFFIX + params.FILE_EXTENSION)).exists()) {
				count++;
				String name = OutputSettings.PATHWAYS_FILE;
				if (count == 1) {
					//                    int index = name.lastIndexOf(".");
					//                   KPMParameters.PATHWAYS_FILE = name.substring(0, index) + "(" + String.valueOf(count) + ")"
					//                            + name.substring(index, name.length());
					OutputSettings.PATHWAYS_FILE += "(" + String.valueOf(count) + ")";
				} else {
					int aux = count - 1;
					OutputSettings.PATHWAYS_FILE = name.replace("(" + aux + ")", "(" + count + ")");
				}

			}
		}

		if (OutputSettings.PATHWAYS_STATS_FILE != null) {
			int count = 0;
			while ((new File(OutputSettings.PATHWAYS_STATS_FILE + "-"
					+ params.SUFFIX + params.FILE_EXTENSION)).exists()) {
				count++;
				String name = OutputSettings.PATHWAYS_STATS_FILE;
				if (count == 1) {
					OutputSettings.PATHWAYS_STATS_FILE += "(" + String.valueOf(count) + ")";
				} else {
					int aux = count - 1;
					OutputSettings.PATHWAYS_STATS_FILE = name.replace("(" + aux + ")", "(" + count + ")");
				}
			}
		}

		if (OutputSettings.GENE_STATS_FILE != null) {
			int count = 0;
			while ((new File(OutputSettings.GENE_STATS_FILE + "-"
					+ params.SUFFIX + params.FILE_EXTENSION)).exists()) {
				count++;
				String name = OutputSettings.GENE_STATS_FILE;
				if (count == 1) {
					OutputSettings.GENE_STATS_FILE += "(" + String.valueOf(count) + ")";
				} else {
					int aux = count - 1;
					OutputSettings.GENE_STATS_FILE = name.replace("(" + aux + ")", "(" + count + ")");
				}

			}
		}

		if (OutputSettings.DATASETS_STATS_FILE != null) {
			int count = 0;
			while ((new File(OutputSettings.DATASETS_STATS_FILE + "-"
					+ params.SUFFIX + params.FILE_EXTENSION)).exists()) {
				count++;
				String name = OutputSettings.DATASETS_STATS_FILE;
				if (count == 1) {
					OutputSettings.DATASETS_STATS_FILE += "(" + String.valueOf(count) + ")";
				} else {
					int aux = count - 1;
					OutputSettings.DATASETS_STATS_FILE = name.replace("(" + aux + ")", "(" + count + ")");
				}

			}
		}

		return params;
	}
}
