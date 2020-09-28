package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.mpg.mpiinf.ag1.kpm.main.Main;

public class PriorityFileParser {
	public static Set<String> parse(String file) {
		BufferedReader br = null;
		Set<String> ret = new HashSet<String>();
		try {
			if (!new File(file).isFile()) {
				return ret;
			}
			br = new BufferedReader(new FileReader(file));
			String line;
			while ((line = br.readLine()) != null) {
				String[] entries = line.split("\t");
				String id = entries[0].trim();
				ret.add(id);
			}
			br.close();
			return ret;
		} catch (IOException ioe) {
			Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ioe);
		} finally {
			try {
				br.close();
			} catch (IOException ex) {
				Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
			}
		}

		return ret;
	}
}
