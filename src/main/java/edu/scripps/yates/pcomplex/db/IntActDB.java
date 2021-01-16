package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import gnu.trove.map.hash.THashMap;

public class IntActDB extends ProteinComplexDB {
	private final static Logger log = Logger.getLogger(IntActDB.class);

	/**
	 * 
	 * @param inputFile
	 * @param organismTaxID
	 * @param uplr
	 * @throws IOException
	 */
	public IntActDB(File inputFile, String organismTaxID, UniprotProteinLocalRetriever uplr) throws IOException {

		super(null, "IntAct - " + organismTaxID, false, uplr);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String line = null;
			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			int numLine = 0;
			int numInOrganism = 0;
			log.info("Reading IntAct database for " + organismTaxID + " from file " + inputFile.getAbsolutePath());
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
					final String taxID = split[indexesByHeaders.get("Taxonomy identifier")];

					if (organismTaxID == null || taxID.equals(organismTaxID)) {
						numInOrganism++;
						final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID);
						// IMPORTANT
						// because some complexes will have the same components than others and ALL will
						// have an id, we use that ID to differentiate them in the Sets (affecting
						// hashcode and equals methods)
						proteinComplex.setUseIdInEqualsMethod(true);
						proteinComplex.setName(proteinComplexName);
						proteinComplex.setOrganism(taxID);
						final String componentsAccString = split[indexesByHeaders
								.get("Identifiers (and stoichiometry) of molecules in complex")];
						// it could be like P84022(1)|Q13485(1)|Q15796(1) or
						// O14746(2)|URS00004A7003_9606(2)
						// CHEBI:29035(1)|CHEBI:29105(0)|CHEBI:29108(0)|P05109(2)|P06702(2)
						final String[] split2 = componentsAccString.split("\\|");
						for (int i = 0; i < split2.length; i++) {
							final String rawACC = split2[i];
							// it could be like CHEBI:29035(1) or O14746(2) or URS00004A7003_9606(2)
							// first we remove the stoichiometry
							final String regexp = "(\\S+)(\\(\\d+\\)$)";
							final Pattern pattern = Pattern.compile(regexp);
							final Matcher matcher = pattern.matcher(rawACC);
							String acc = null;
							if (matcher.find()) {
								acc = matcher.group(1);

							} else {
								acc = rawACC;
							}
							final ProteinComponent component = new ProteinComponent(acc, null);
							proteinComplex.addComponent(component);
						}
						proteinComplexes.add(proteinComplex);

					} else {
						System.out.println(taxID);
					}
				}
			}
			setReady();
			log.info("IntAct loaded for organism " + organismTaxID + ": " + getProteinComplexes().size()
					+ " complexes read in " + numLine + " lines");
		} finally {
			br.close();
		}
	}
}
