package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComplexExistenceCriteria;
import edu.scripps.yates.pcomplex.util.PComplexUtil;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinComplexDB {
	protected final List<ProteinComplex> proteinComplexes = new ArrayList<ProteinComplex>();
	private final Map<String, Set<ProteinComplex>> proteinComplexesByComponent = new THashMap<String, Set<ProteinComplex>>();
	private final String name;
	private final static Logger log = Logger.getLogger(ProteinComplexDB.class);

	public ProteinComplexDB(File inputFile, String name, boolean load) throws IOException {
		this(inputFile, name, null, load);
	}

	public ProteinComplexDB(File inputFile, String name, File mappingFile, boolean load) throws IOException {
		this.name = name;
		if (load) {
			Map<String, String> geneToUniprotMap = null;
			if (mappingFile != null) {
				geneToUniprotMap = new THashMap<String, String>();
				BufferedReader br2 = null;
				try {
					br2 = new BufferedReader(new FileReader(mappingFile));
					String line = null;
					int numline = 0;
					while ((line = br2.readLine()) != null) {
						numline++;
						if (numline == 1) {
							continue;
						}
						if (line.startsWith("\"")) {
							continue;
						}
						final String[] split = line.split(",");
						final String geneName = split[3];
						final String uniprotACC = split[0];
						geneToUniprotMap.put(geneName, uniprotACC);
					}

				} finally {
					br2.close();
				}
			}
			BufferedReader br = null;
			try {

				br = new BufferedReader(new FileReader(inputFile));
				String line = null;
				int numLine = 1;
				while ((line = br.readLine()) != null) {
					final ProteinComplex proteinComplex = new ProteinComplex(String.valueOf(numLine++), true);
					proteinComplexes.add(proteinComplex);
					final String[] split = line.split("\t");
					for (final String geneName : split) {
						if (geneToUniprotMap != null && !geneToUniprotMap.containsKey(geneName)) {
							log.warn("Mapping for gene name: " + geneName + " is not found");
							addComponentToProteinComplex(proteinComplex, geneName);
							continue;
						}
						if (geneToUniprotMap != null) {
							addComponentToProteinComplex(proteinComplex, geneToUniprotMap.get(geneName));
						} else {
							addComponentToProteinComplex(proteinComplex, geneName);
						}
					}
				}
				setReady();
			} finally {
				br.close();
			}

		}
	}

	protected void addComponentToProteinComplex(ProteinComplex proteinComplex, String uniprotACC) {
		proteinComplex.addComponent(uniprotACC);

	}

	protected void setReady() {

		for (final ProteinComplex proteinComplex : proteinComplexes) {
			for (final String prot : proteinComplex.getComponents()) {
				if (proteinComplexesByComponent.containsKey(prot)) {
					proteinComplexesByComponent.get(prot).add(proteinComplex);
				} else {
					final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
					set.add(proteinComplex);
					proteinComplexesByComponent.put(prot, set);
				}
			}
		}

	}

	public List<ProteinComplex> getProteinComplexes() {
		return proteinComplexes;
	}

	public Set<ProteinComplex> getProteinComplexesByProtein(Protein protein) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		if (proteinComplexesByComponent.containsKey(protein.getAcc())) {
			ret.addAll(proteinComplexesByComponent.get(protein.getAcc()));
		}
		if (proteinComplexesByComponent.containsKey(protein.getGene())) {
			ret.addAll(proteinComplexesByComponent.get(protein.getGene()));
		}
		return ret;
	}

	public Set<ProteinComplex> getProteinComplexesByProtein(String proteinKey) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		if (proteinComplexesByComponent.containsKey(proteinKey)) {
			ret.addAll(proteinComplexesByComponent.get(proteinKey));
		}

		return ret;
	}

	/**
	 * Get the set of {@link ProteinComplex}es that are fully contained in the input
	 * protein collection
	 * 
	 * @param proteinsAndGenes
	 * @param minNumComponents minimum number of components necessary to be counted
	 * @return
	 */
	public Set<ProteinComplex> getCompletedProteinComplexes(Set<String> proteinsAndGenes) {
		return getCompleteProteinComplexes(proteinsAndGenes, 0);
	}

	/**
	 * Get the set of {@link ProteinComplex}es that are fully contained in the input
	 * protein collection
	 * 
	 * @param proteinsAndGenes
	 * @param minNumComponents minimum number of components necessary to be counted
	 * @return
	 */
	public Set<ProteinComplex> getCompleteProteinComplexes(Set<String> proteinsAndGenes, int minNumComponents) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();

		for (final ProteinComplex proteinComplex : proteinComplexes) {
			if (proteinComplex.getComponents().size() >= minNumComponents) {
				if (proteinComplex.isFullyRepresented(proteinsAndGenes)) {
					ret.add(proteinComplex);
				}
			}
		}

		return ret;
	}

	/**
	 * Gets the {@link ProteinComplex} list from this {@link ProteinComplexDB}
	 * following a {@link ProteinComplexExistenceCriteria}
	 * 
	 * @param proteinsAndGenes
	 * @param existenceCriteria
	 * @return
	 */
	public Set<ProteinComplex> getCompleteProteinComplexes(Set<String> proteinsAndGenes,
			ProteinComplexExistenceCriteria existenceCriteria, int minNumComponents) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();

		for (final ProteinComplex proteinComplex : proteinComplexes) {
			if (proteinComplex.getComponents().size() >= minNumComponents) {
				if (existenceCriteria.considersExisting(proteinComplex, proteinsAndGenes)) {

					ret.add(proteinComplex);

				}
			}
		}

		return ret;
	}

	/**
	 * Get the set of {@link ProteinComplex}es that are fully contained in the input
	 * protein collection
	 * 
	 * @param proteins
	 * @return
	 */
	public Set<ProteinComplex> getPartiallyCompletedProteinComplexes(Collection<Protein> proteins) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		final Set<String> keys = PComplexUtil.getKeysFromProteins(proteins);

		for (final ProteinComplex proteinComplex : proteinComplexes) {
			if (proteinComplex.isPartiallyRepresented(keys) && !proteinComplex.isFullyRepresented(keys)) {
				ret.add(proteinComplex);
			}
		}

		return ret;
	}

	public String getName() {
		return name;
	}

	/**
	 * Returns the set of proteins that are found to be part of protein complexes
	 * that are completed among the proteins
	 * 
	 * @param proteins
	 * @param minNumComponentsInComplex
	 * @return
	 */
	public Set<Protein> getProteinsInCompleteProteinComplexes(Set<Protein> proteins, int minNumComponentsInComplex) {
		final Set<Protein> ret = new THashSet<Protein>();
		// get the protein complexes that are completed represented in the input
		// proteins
		final Set<String> proteinsAndGenes = new THashSet<String>();
		for (final Protein protein : proteins) {
			proteinsAndGenes.add(protein.getAcc());
			if (protein.getGene() != null) {
				proteinsAndGenes.add(protein.getGene());
			}
		}
		final Set<ProteinComplex> proteinComplexes = getCompleteProteinComplexes(proteinsAndGenes,
				minNumComponentsInComplex);
		// get the map that gets the proteins per key names of them
		final Map<String, Set<Protein>> keys = PComplexUtil.getKeyMapFromProteins(proteins);
		// get the key names of the protein complexes
		final Set<String> proteinComplexesComponents = PComplexUtil.getKeysFromProteinComplexes(proteinComplexes);
		// look all the key names of the proteins and see if they are in the
		// protein complexes
		for (final String proteinKey : keys.keySet()) {
			if (proteinComplexesComponents.contains(proteinKey)) {
				ret.addAll(keys.get(proteinKey));
			}
		}
		return ret;
	}

	/**
	 * Returns the set of proteins that are found to be part of protein complexes,
	 * not necessarily covering all their members of the complex (partial complexes)
	 * 
	 * @param proteins
	 * @return
	 */
	public Set<Protein> getProteinsInPartiallyCompletedProteinComplexes(Set<Protein> proteins) {
		final Set<Protein> ret = new THashSet<Protein>();
		// get the protein complexes from which we have at least one component
		// in the input proteins
		final Set<ProteinComplex> proteinComplexes = getPartiallyCompletedProteinComplexes(proteins);
		// get the map that gets the proteins per key names of them
		final Map<String, Set<Protein>> keys = PComplexUtil.getKeyMapFromProteins(proteins);
		// get the key names of the protein complexes
		final Set<String> proteinComplexesComponents = PComplexUtil.getKeysFromProteinComplexes(proteinComplexes);
		// look all the key names of the proteins and see if they are in the
		// protein complexes
		for (final String proteinKey : keys.keySet()) {
			if (proteinComplexesComponents.contains(proteinKey)) {
				ret.addAll(keys.get(proteinKey));
			}
		}
		return ret;

	}

	public Set<Protein> getProteinsInCompleteProteinComplexes(Set<Protein> proteins,
			ProteinComplexExistenceCriteria existenceCriteria, int minNumComponents) {
		final Set<Protein> ret = new THashSet<Protein>();
		// get the protein complexes that are completed represented in the input
		// proteins
		final Set<String> proteinsAndGenes = new THashSet<String>();
		for (final Protein protein : proteins) {
			proteinsAndGenes.add(protein.getAcc());
			if (protein.getGene() != null) {
				proteinsAndGenes.add(protein.getGene());
			}
		}
		final Set<ProteinComplex> proteinComplexes = getCompleteProteinComplexes(proteinsAndGenes, existenceCriteria,
				minNumComponents);
		// get the map that gets the proteins per key names of them
		final Map<String, Set<Protein>> keys = PComplexUtil.getKeyMapFromProteins(proteins);
		// get the key names of the protein complexes
		final Set<String> proteinComplexesComponents = PComplexUtil.getKeysFromProteinComplexes(proteinComplexes);
		// look all the key names of the proteins and see if they are in the
		// protein complexes
		for (final String proteinKey : keys.keySet()) {
			if (proteinComplexesComponents.contains(proteinKey)) {
				ret.addAll(keys.get(proteinKey));
			}
		}
		return ret;
	}
}
