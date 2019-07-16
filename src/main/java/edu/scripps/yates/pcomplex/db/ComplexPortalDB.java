package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import gnu.trove.map.hash.THashMap;

public class ComplexPortalDB extends ProteinComplexDB {

	public ComplexPortalDB(File inputFile) throws IOException {
		super(null, "Complex Portal", false);
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
					final String proteinComplexID = split[indexesByHeaders.get("#Complex ac")];
					final String proteinComplexName = split[indexesByHeaders.get("Recommended name")];
					final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID, false);
					proteinComplexes.add(proteinComplex);
					proteinComplex.setName(proteinComplexName);
					final String componentsString = split[indexesByHeaders
							.get("Identifiers (and stoichiometry) of molecules in complex")];
					final String[] split2 = componentsString.split("\\|");
					for (String component : split2) {
						// component is like Q07812(8)
						component = component.substring(0, component.indexOf("("));
						addComponentToProteinComplex(proteinComplex, component);
					}
				}
			}
			setReady();
		} finally {
			br.close();
		}
	}

}
