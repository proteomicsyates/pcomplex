package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import gnu.trove.map.hash.THashMap;

public class CoreCorumDB extends ProteinComplexDB {

	public CoreCorumDB(File inputFile, String organism) throws IOException {
		super(null, "Core complexes - CORUM", false);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String line = null;
			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			int numLine = 0;
			while ((line = br.readLine()) != null) {
				numLine++;
				final String[] split = line.split("\t");
				if (numLine == 1) {
					for (int index = 0; index < split.length; index++) {
						indexesByHeaders.put(split[index], index);
					}
				} else {
					final String proteinComplexID = split[indexesByHeaders.get("ComplexID")];
					final String proteinComplexName = split[indexesByHeaders.get("ComplexName")];
					final String organism2 = split[indexesByHeaders.get("Organism")];
					if (organism == null || organism2.equals(organism)) {
						final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID, false);
						proteinComplexes.add(proteinComplex);
						proteinComplex.setName(proteinComplexName);
						final String componentsString = split[indexesByHeaders.get("subunits(UniProt IDs)")];
						final String[] split2 = componentsString.split(";");
						for (final String component : split2) {
							addComponentToProteinComplex(proteinComplex, component);
						}
					}
				}
			}
			setReady();
		} finally {
			br.close();
		}
	}

}
