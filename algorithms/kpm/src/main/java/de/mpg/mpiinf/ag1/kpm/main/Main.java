/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.main;

import de.mpg.mpiinf.ag1.kpm.KPMRunHandler;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.Program;
import de.mpg.mpiinf.ag1.kpm.parsers.ArgsParametersParser;
import de.mpg.mpiinf.ag1.kpm.parsers.InputFileParser;
import de.mpg.mpiinf.ag1.kpm.parsers.PropertiesFileParser;
import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathAlgorithms;
import dk.sdu.kpm.KPMSettings;

/**
 * @author nalcaraz
 */

public class Main {

    public static void main(String[] args) {
        Main main = new Main();
        main.run(args);
    }

    //Runner method for calling KeyPathwayMiner
    public void run(String[] args) {
        String datasetFolder = "resources/";
        String propertiesFile = "kpm.properties";

        System.out.println("Running KeyPathWayMiner standalone.");
        System.out.println("Properties file:" + propertiesFile);
        System.out.println("Default datasets from: " + datasetFolder);
        kpm(args, datasetFolder, propertiesFile);
    }

    //Runner method for calling KeyPathwayMiner from R
    public void runR(String[] args, String datasetFolder, String propertiesFile) {
        System.out.println("Running KeyPathWayMiner standalone via R.");
        System.out.println("Properties file: " + propertiesFile);
        System.out.println("Default datasets from: " + datasetFolder);
        kpm(args, datasetFolder, propertiesFile);

    }

    // Prepares objects, parses arguments and runs KeyPathwayMiner
    private void kpm(String[] args, String datasetFolder, String propertiesFile) {
        //KPM default settings
        KPMSettings kpmSettings = new KPMSettings();

        //Parse parameters from the kpm.properties file, initialize default parameters and runs KeyPathwayMiner
        PropertiesFileParser pfp = new PropertiesFileParser(kpmSettings, propertiesFile);

        Parameters params = pfp.parse(datasetFolder);

        try {
            //Read commandline arguments
            ArgsParametersParser argsParser = new ArgsParametersParser(kpmSettings);
            params = argsParser.parse(args, params);

            //Check if all input files are correct
            InputFileParser ifp = new InputFileParser(kpmSettings);
            params = ifp.parse(params);

            //Start KeyPathwayMiner
            if (params.PROGRAM == Program.SP) {
                System.out.println("PROGRAM SELECTED: Shortest Paths");
                ShortestPathAlgorithms.shortestPathways(params.GRAPH_FILE, params);
            } else {
                KPMRunHandler kpmHandler = new KPMRunHandler(kpmSettings);
                kpmHandler.runBatch(params);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


}
