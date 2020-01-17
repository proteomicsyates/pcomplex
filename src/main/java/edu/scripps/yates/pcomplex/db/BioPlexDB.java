package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class BioPlexDB extends ProteinComplexDB {
	private final static Logger log = Logger.getLogger(BioPlexDB.class);

	public BioPlexDB(List<File> inputFiles, UniprotProteinLocalRetriever uplr) throws IOException {
		super(null, "BioPlex v1-3", false, uplr);
		for (final File inputFile : inputFiles) {
			log.info("Reading file " + inputFile.getAbsolutePath());
			BufferedReader br = null;
			String line = null;
			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			try {
				br = new BufferedReader(new FileReader(inputFile));
				indexesByHeaders.clear();

				int numLine = 0;
				while ((line = br.readLine()) != null) {
					numLine++;
					final String[] split = line.split("\t");
					if (numLine == 1) {
						for (int index = 0; index < split.length; index++) {
							indexesByHeaders.put(split[index].replace(" ", ""), index);
						}
					} else {
						final int geneAID = Integer.valueOf(split[indexesByHeaders.get("GeneA")]);
						final int geneBID = Integer.valueOf(split[indexesByHeaders.get("GeneB")]);
						// returns the concatenated String of the two numbers after sorting them
						final String proteinComplexID = getComplexID(geneAID, geneBID);
						final String proteinComplexName = proteinComplexID;
						final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID);

						proteinComplex.setName(proteinComplexName);

						// component A
						final String accA = split[indexesByHeaders.get("UniprotA")];
						final String geneA = split[indexesByHeaders.get("SymbolA")];
						final ProteinComponent pcA = new ProteinComponent(accA, geneA);
						proteinComplex.addComponent(pcA);

						// component B
						final String accB = split[indexesByHeaders.get("UniprotB")];
						final String geneB = split[indexesByHeaders.get("SymbolB")];
						final ProteinComponent pcB = new ProteinComponent(accB, geneB);
						proteinComplex.addComponent(pcB);

						proteinComplexes.add(proteinComplex);
					}
				}

			} catch (final Exception e) {
				e.printStackTrace();
				log.error("Error in line : " + line);
			} finally {
				br.close();
			}
		}
		// look for uniprot accessions using the gene, when they are null
		final Set<String> accs = new THashSet<String>();
		for (final ProteinComplex complex : proteinComplexes) {
			final List<ProteinComponent> componentList = complex.getComponentList();
			for (final ProteinComponent component : componentList) {
				if ("UNKNOWN".equals(component.getAcc())) {
					final Set<String> uniprots = ProteinComplexAnalyzer.getUniprotGeneMapping()
							.mapGeneToUniprotACC(component.getGene());
					if (uniprots != null && !uniprots.isEmpty()) {
						if (uniprots.size() == 1) {
							component.setAcc(uniprots.iterator().next());
						} else {
							accs.addAll(uniprots);
						}
					}
				}
			}
		}
		// get all accs first
		final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
		for (final ProteinComplex complex : proteinComplexes) {
			final List<ProteinComponent> componentList = complex.getComponentList();
			for (final ProteinComponent component : componentList) {
				if ("UNKNOWN".equals(component.getAcc())) {
					final Set<String> uniprots = ProteinComplexAnalyzer.getUniprotGeneMapping()
							.mapGeneToUniprotACC(component.getGene());
					if (uniprots.isEmpty()) {
						continue;
					}
					final String acc = getReviewedIfPossible(uniprots, uplr);
					component.setAcc(acc);
				}
			}
		}
		setReady();
	}

	private String getReviewedIfPossible(Set<String> uniprots, UniprotProteinLocalRetriever uplr) {
		final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, uniprots);
		for (final String acc : uniprots) {
			if (annotatedProteins.containsKey(acc)) {
				final Entry entry = annotatedProteins.get(acc);
				if (UniprotEntryUtil.isSwissProt(entry)) {
					return acc;
				}
			}
		}
		return uniprots.iterator().next();
	}

	private String getComplexID(int geneAID, int geneBID) {
		if (geneAID < geneBID) {
			return "" + geneAID + "" + geneBID;
		} else {
			return "" + geneBID + "" + geneAID;
		}
	}

}
