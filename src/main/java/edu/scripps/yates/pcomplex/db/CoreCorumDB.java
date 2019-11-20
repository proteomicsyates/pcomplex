package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.mi.MolecularInteractionsOntologyClient;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import gnu.trove.map.hash.THashMap;

public class CoreCorumDB extends ProteinComplexDB {
	private final static Logger log = Logger.getLogger(CoreCorumDB.class);
	private static final Map<String, Boolean> isValid = new THashMap<String, Boolean>();
	private final static String validityFileName = "valid.txt";
	private boolean isloaded = false;
	private final File inputFile;

	public CoreCorumDB(File inputFile, String organism, UniprotProteinLocalRetriever uplr,
			boolean filterByPurificationMethod) throws IOException {

		super(null, "Core complexes - CORUM", false, uplr);
		this.inputFile = inputFile;
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String line = null;
			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			int numLine = 0;
			int numInOrganism = 0;
			int numValid = 00;
			log.info("Reading CORUM database for " + organism + " from file " + inputFile.getAbsolutePath());
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
					final String purificationMethod = split[indexesByHeaders
							.get("Protein complex purification method")];

					if (organism == null || organism2.equals(organism)) {
						numInOrganism++;
						if (!filterByPurificationMethod || isValid(purificationMethod)) {
							numValid++;
							final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID);

							proteinComplex.setName(proteinComplexName);
							proteinComplex.setOrganism(organism2);
							final String componentsAccString = split[indexesByHeaders.get("subunits(UniProt IDs)")];
							final String[] split2 = componentsAccString.split(";");
							final String componentsGeneString = split[indexesByHeaders.get("subunits(Gene name)")];
							final String[] split3 = componentsGeneString.split(";");
							final String componentsNamesString = split[indexesByHeaders.get("subunits(Protein name)")];
							final String[] split4 = componentsNamesString.split(";");
							for (int i = 0; i < split2.length; i++) {
								final String acc = split2[i];
								String gene = null;
								if (split3.length > i) {
									gene = split3[i];
								} else {
									log.info(componentsGeneString);
								}
								String proteinName = null;
								if (split4.length > i) {
									proteinName = split4[i];
								} else {
									log.info(componentsGeneString);
								}

								final ProteinComponent component = new ProteinComponent(acc, gene);
								component.setProteinName(proteinName);
								proteinComplex.addComponent(component);
							}
							proteinComplexes.add(proteinComplex);
						}
					}
				}
			}
			setReady();
			log.info("CORUM loaded for organism " + organism);
			log.info(numValid + " complexes are valid as biochemical ones");
			log.info((numInOrganism - numValid) + " complexes discarded as non biochemical ones");
		} finally {
			br.close();
		}
	}

	private boolean isValid(String purificationMethod) throws IOException {
		if (!isloaded) {
			loadIsValid();
		}
		if (isValid.containsKey(purificationMethod)) {
			return isValid.get(purificationMethod);
		}
		final List<String> termIDs = getTermIDsFromPurificationMethod(purificationMethod);
		for (final String termID : termIDs) {
			if (termID == null) {
				continue;
			}
			if (isValid.containsKey(termID)) {
				return isValid.get(termID);
			}
			final boolean isBiochemical = MolecularInteractionsOntologyClient.containsBiochemicalAsParent(termID);
			addToValidity(termID, isBiochemical);

			if (isBiochemical) {
				isValid.put(purificationMethod, true);
				return true;
			}
		}

		return false;
	}

	private void addToValidity(String termID, boolean isValidBoolean) throws IOException {
		isValid.put(termID, isValidBoolean);
		// now append to file
		final FileWriter fw = new FileWriter(getIsValidFile(), true);
		fw.write(termID + "\t" + isValidBoolean + "\n");
		fw.close();
	}

	private void loadIsValid() throws IOException {
		final File isValidFile = getIsValidFile();
		if (isValidFile.exists()) {
			final List<String> lines = Files.readAllLines(isValidFile.toPath());
			for (final String line : lines) {
				final String[] split = line.split("\t");
				final String term = split[0];
				final boolean valid = Boolean.valueOf(split[1]);
				isValid.put(term, valid);
			}
			log.info(isValid.size() + " validities loaded from file " + isValidFile.getAbsolutePath());
		}
		isloaded = true;
	}

	private File getIsValidFile() {
		return new File(inputFile.getParentFile().getAbsolutePath() + File.separator + validityFileName);
	}

	private List<String> getTermIDsFromPurificationMethod(String purificationMethod) {
		final List<String> purificationMethods = new ArrayList<String>();
		if (purificationMethod.contains(";")) {
			final String[] split = purificationMethod.split(";");
			for (final String string : split) {
				purificationMethods.add(string);
			}
		} else {
			purificationMethods.add(purificationMethod);
		}
		return purificationMethods.stream().map(pm -> getID(pm)).collect(Collectors.toList());
	}

	private String getID(String pm) {
		if (pm.contains("-")) {
			final String substring = pm.substring(0, pm.indexOf("-"));
			return substring;
		} else {
			return null;
		}
	}
}
