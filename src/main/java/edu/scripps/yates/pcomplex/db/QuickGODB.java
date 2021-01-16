package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import gnu.trove.map.hash.THashMap;

public class QuickGODB extends ProteinComplexDB {
	private final static Logger log = Logger.getLogger(QuickGODB.class);

	/**
	 * 
	 * @param inputFile
	 * @param organismTaxID
	 * @param uplr
	 * @throws IOException
	 */
	public QuickGODB(File inputFile, String organismTaxID, UniprotProteinLocalRetriever uplr) throws IOException {

		super(null, "QuickGO - " + organismTaxID, false, uplr);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String line = null;
			// in this case, we have a table in which in each line is a component of a
			// complex, so that we have to keep the complex in a map by its ID so that we
			// can add its elements from different rows
			// we have to have this map here because the one in the super class is only
			// ready when the file has read entirely

			final Map<String, ProteinComplex> complexesByID = new THashMap<String, ProteinComplex>();

			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			int numLine = 0;
			int numInOrganism = 0;
			log.info("Reading QuickGO database for " + organismTaxID + " from file " + inputFile.getAbsolutePath());
			while ((line = br.readLine()) != null) {
				numLine++;
				final String[] split = line.split("\t");
				if (numLine == 1) {
					for (int index = 0; index < split.length; index++) {
						indexesByHeaders.put(split[index], index);
					}
				} else {
					final String proteinComplexID = split[indexesByHeaders.get("GO TERM")];
					final String taxID = split[indexesByHeaders.get("TAXON ID")];

					if (organismTaxID == null || taxID.equals(organismTaxID)) {
						numInOrganism++;
						ProteinComplex proteinComplex = null;
						if (complexesByID.containsKey(proteinComplexID)) {
							proteinComplex = complexesByID.get(proteinComplexID);
						} else {
							proteinComplex = new ProteinComplex(proteinComplexID);
							// IMPORTANT
							// because some complexes will have the same components than others and ALL will
							// have an id, we use that ID to differentiate them in the Sets (affecting
							// hashcode and equals methods)
							proteinComplex.setUseIdInEqualsMethod(true);
							complexesByID.put(proteinComplexID, proteinComplex);
							proteinComplex.setOrganism(taxID);
						}

						final String acc = split[indexesByHeaders.get("GENE PRODUCT ID")];
						final String gene = split[indexesByHeaders.get("SYMBOL")];

						// usually are uniprot single Accs

//						if (FastaParser.isUniProtACC(acc)) {
						final ProteinComponent component = new ProteinComponent(acc, gene);
						proteinComplex.addComponent(component);
//						} else {
//
//						}

					} else {
						log.warn(line + " skipped because doesn't belongs to taxid " + organismTaxID);
					}
				}
			}
			// add all complexes in the map to the super class
			for (final String complexID : complexesByID.keySet()) {
				final ProteinComplex proteinComplex = complexesByID.get(complexID);
				super.proteinComplexes.add(proteinComplex);
			}
			setReady();
			log.info("QuickGO loaded for organism " + organismTaxID + ": " + complexesByID.size()
					+ " complexes read in " + numLine + " lines");
		} finally {
			br.close();
		}
	}
}
