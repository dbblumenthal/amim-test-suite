package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import au.com.bytecode.opencsv.CSVReader;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.main.Main;
import dk.sdu.kpm.KPMSettings;

public class DatasetsFileParser {

	private volatile KPMSettings kpmSettings;

	public DatasetsFileParser(KPMSettings settings){
		this.kpmSettings = settings;
	}

	public Parameters parse(char separator, boolean hasHeader, Parameters params) {
		try {

            HashMap<String, String> id2path = new HashMap<String, String>();
            HashMap<String, Integer> id2param = new HashMap<String, Integer>();
			CSVReader cr =
					new CSVReader(new FileReader(params.DATASETS_FILE), separator);
			String[] nextLine = cr.readNext();

			int cols = nextLine.length;
			int row = 1;
			if (!hasHeader) {
				String id = "";
				String filePath = "";
				int lParam = 0;
				if (cols == 2) {
					id = String.valueOf(row);
					lParam = Integer.parseInt(nextLine[0]);
					filePath = nextLine[1];
				} else if (cols == 3) {
					id = nextLine[0];
					lParam = Integer.parseInt(nextLine[1]);
					filePath = nextLine[2];
				}
				if (!new File(filePath).isFile()) {
					System.out.println("Matrix file " + filePath + " does not exist !");
					System.exit(-1);
				}
				id2param.put(id, lParam);
				id2path.put(id, filePath);
				row++;
			}
			while ((nextLine = cr.readNext()) != null) {
				String id;
				String filePath;
				int lParam;
				if (nextLine.length < 2) {
					break;
				}
				if (nextLine.length == 2) {
					id = String.valueOf(row);
					lParam = Integer.parseInt(nextLine[0]);
					filePath = nextLine[1];
				} else {
					id = nextLine[0];
					lParam = Integer.parseInt(nextLine[1]);
					filePath = nextLine[2];
				}
				if (!new File(filePath).isFile()) {
					System.out.println("Matrix file " + filePath + " does not exist !");
					System.exit(-1);
				}
				id2param.put(id, lParam);
				id2path.put(id, filePath);

				kpmSettings.MIN_L.put(id, lParam);
				kpmSettings.INC_L.put(id, lParam);
				kpmSettings.MAX_L.put(id, lParam);

				row++;
			}

			kpmSettings.CASE_EXCEPTIONS_MAP = id2param;
			kpmSettings.MATRIX_FILES_MAP = id2path;
			cr.close();
		} catch (IOException ex) {
			Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
		}

		return params;
	}
}
