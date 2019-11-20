package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComplexExistenceCriteria;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import edu.scripps.yates.pcomplex.util.PComplexUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinComplexDB {
	protected final Set<ProteinComplex> proteinComplexes = new THashSet<ProteinComplex>();
	protected final Map<String, Set<ProteinComplex>> proteinComplexesByComponent = new THashMap<String, Set<ProteinComplex>>();
	private final THashMap<String, ProteinComplex> proteinComplexesByID = new THashMap<String, ProteinComplex>();

	private final String name;
	private final static Logger log = Logger.getLogger(ProteinComplexDB.class);
	private final UniprotProteinLocalRetriever uplr;
	private boolean indexByGene;

	public ProteinComplexDB(File inputFile, String name, boolean load, UniprotProteinLocalRetriever uplr)
			throws IOException {
		this(inputFile, name, null, load, uplr);
	}

	public ProteinComplexDB(File inputFile, String name, boolean load, UniprotProteinLocalRetriever uplr,
			boolean indexByGene) throws IOException {
		this(inputFile, name, null, load, uplr, indexByGene);
	}

	public ProteinComplexDB(File inputFile, String name, File mappingFile, boolean load,
			UniprotProteinLocalRetriever uplr) throws IOException {
		this(inputFile, name, mappingFile, load, uplr, true);
	}

	public ProteinComplexDB(File inputFile, String name, File mappingFile, boolean load,
			UniprotProteinLocalRetriever uplr, boolean indexByGene) throws IOException {
		this.uplr = uplr;
		this.name = name;
		this.indexByGene = indexByGene;
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
					final ProteinComplex proteinComplex = new ProteinComplex(String.valueOf(numLine++));

					final String[] split = line.split("\t");
					for (final String geneName : split) {
						if (geneToUniprotMap != null && !geneToUniprotMap.containsKey(geneName)) {
							log.warn("Mapping for gene name: " + geneName + " is not found");
							final ProteinComponent pc = new ProteinComponent(null, geneName);
							proteinComplex.addComponent(pc);
							continue;
						}
						if (geneToUniprotMap != null) {
							final String acc = geneToUniprotMap.get(geneName);
							final ProteinComponent pc = new ProteinComponent(acc, geneName);
							proteinComplex.addComponent(pc);
						} else {
							final ProteinComponent pc = new ProteinComponent(null, geneName);
							proteinComplex.addComponent(pc);
						}
					}
					proteinComplexes.add(proteinComplex);
				}
				setReady();
			} finally {
				br.close();
			}

		}
	}

	protected void setReady() {
		proteinComplexesByID.clear();
		proteinComplexesByComponent.clear();
		// get all annotations and set protein names
		final Set<String> accs = new THashSet<String>();
		for (final ProteinComplex proteinComplex : proteinComplexes) {
			proteinComplexesByID.put(proteinComplex.getId(), proteinComplex);
			for (final ProteinComponent prot : proteinComplex.getComponentList()) {
				if (prot.getAcc() != null) {
					accs.add(prot.getAcc());
					if (proteinComplexesByComponent.containsKey(prot.getAcc())) {
						proteinComplexesByComponent.get(prot.getAcc()).add(proteinComplex);
					} else {
						final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
						set.add(proteinComplex);
						proteinComplexesByComponent.put(prot.getAcc(), set);
					}
				}
				if (indexByGene && prot.getGene() != null) {
					if (proteinComplexesByComponent.containsKey(prot.getGene())) {
						proteinComplexesByComponent.get(prot.getGene()).add(proteinComplex);
					} else {
						final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
						set.add(proteinComplex);
						proteinComplexesByComponent.put(prot.getGene(), set);
					}
				}
			}
		}
		if (!accs.isEmpty()) {
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			for (final ProteinComplex proteinComplex : proteinComplexes) {
				for (final ProteinComponent prot : proteinComplex.getComponentList()) {
					if (prot.getAcc() != null && prot.getProteinName() == null) {
						if (annotatedProteins.containsKey(prot.getAcc())) {
							final Entry entry = annotatedProteins.get(prot.getAcc());
							final String proteinDescription = UniprotEntryUtil.getProteinDescription(entry);
							prot.setProteinName(proteinDescription);
						}
					}
				}
			}
		}
	}

	public boolean isIndexByGene() {
		return indexByGene;
	}

	public void setIndexByGene(boolean indexByGene) {
		this.indexByGene = indexByGene;
	}

	public Set<ProteinComplex> getProteinComplexes() {
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

	public Set<ProteinComplex> getProteinComplexesByProteins(Collection<String> proteinKeys) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		for (final String proteinKey : proteinKeys) {

			if (proteinComplexesByComponent.containsKey(proteinKey)) {
				ret.addAll(proteinComplexesByComponent.get(proteinKey));
			}
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
			if (proteinComplex.getComponentSet().size() >= minNumComponents) {
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
			if (proteinComplex.getComponentSet().size() >= minNumComponents) {
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

	/**
	 * Merge complexes that have an overlapping score greater than the provided in
	 * the parameter
	 * 
	 * @param maxOverlapScore
	 */
	public void mergeComplexes(double maxOverlapScore) {
		log.info("Merging complexes from REFERENCE set '" + getName() + "' that overlap with overlap score > "
				+ maxOverlapScore);
		final int initialSize = proteinComplexes.size();
		int round = 1;
		while (true) {
			final List<Pair<ProteinComplex, ProteinComplex>> toMerge = new ArrayList<Pair<ProteinComplex, ProteinComplex>>();
			final List<ProteinComplex> proteinComplexList = new ArrayList<ProteinComplex>();
			proteinComplexList.addAll(proteinComplexes);

			for (int i = 0; i < proteinComplexList.size(); i++) {
				final ProteinComplex complex1 = proteinComplexList.get(i);
				for (int j = i + 1; j < proteinComplexList.size(); j++) {
					final ProteinComplex complex2 = proteinComplexList.get(j);
					final double overlap = ClusterEvaluation.getOverlap(complex1, complex2);
					if (overlap > maxOverlapScore) {
						final Pair<ProteinComplex, ProteinComplex> pair = new Pair<ProteinComplex, ProteinComplex>(
								complex1, complex2);
						toMerge.add(pair);
						if (!proteinComplexes.contains(complex1)) {
							System.out.println("WHAT");
						}
						if (!proteinComplexes.contains(complex2)) {
							System.out.println("WHAT");
						}
					}
				}
			}
			log.info(toMerge.size() + " to merge in round " + round);
			// delete te ones in toMerge
			for (final Pair<ProteinComplex, ProteinComplex> pair : toMerge) {
				log.info("Merging:\n" + pair.getFirstelement() + "\n" + pair.getSecondElement());
				boolean removed = proteinComplexes.remove(pair.getFirstelement());
				if (!removed) {
					System.out.println(proteinComplexes.contains(pair.getFirstelement()));
				}
				removed = proteinComplexes.remove(pair.getSecondElement());
				if (!removed) {
					System.out.println(proteinComplexes.contains(pair.getSecondElement()));
				}
			}
			// merge the ones in toMerge
			for (final Pair<ProteinComplex, ProteinComplex> pair : toMerge) {
				final ProteinComplex mergedComplex = new ProteinComplex(pair.getFirstelement().getId());
				mergedComplex.setName(
						"Merged from: " + pair.getFirstelement().getId() + " " + pair.getFirstelement().getName()
								+ " and " + pair.getSecondElement().getId() + " " + pair.getSecondElement().getName());
				for (final ProteinComponent pc : pair.getFirstelement().getComponentList()) {
					mergedComplex.addComponent(pc);
				}
				for (final ProteinComponent pc : pair.getSecondElement().getComponentList()) {
					mergedComplex.addComponent(pc);
				}
				// log.info("\nMerged from:\n" + pair.getFirstelement() + "\n" +
				// pair.getSecondElement() + "\nResult: \t"
				// + mergedComplex);
				proteinComplexes.add(mergedComplex);
			}
			log.info(proteinComplexes.size() + " complexes after round " + round);
			if (toMerge.isEmpty()) {
				break;
			}
			round++;
		}
		proteinComplexesByComponent.clear();
		setReady();
		log.info((initialSize - proteinComplexes.size()) + " complexes merged. Now we have " + proteinComplexes.size()
				+ " complexes (before " + initialSize + ")");
	}

	public void filterByComplexSize(int minComplexSizeInReferenceSetForLearning,
			int maxComplexSizeInReferenceSetForLearning) {
		log.info("Deleting complexes from REFERENCE set '" + getName() + "' with more than "
				+ maxComplexSizeInReferenceSetForLearning + " components and less than "
				+ minComplexSizeInReferenceSetForLearning + " components");
		final int initialSize = proteinComplexes.size();
		int big = 0;
		int small = 0;
		final Iterator<ProteinComplex> iterator = proteinComplexes.iterator();
		while (iterator.hasNext()) {
			final ProteinComplex complex = iterator.next();
			if (complex.size() > maxComplexSizeInReferenceSetForLearning) {
				iterator.remove();
				big++;
			} else if (complex.size() < minComplexSizeInReferenceSetForLearning) {
				iterator.remove();
				small++;
			}
		}
		proteinComplexesByComponent.clear();
		setReady();
		log.info((initialSize - proteinComplexes.size()) + " complexes deleted. Now we have " + proteinComplexes.size()
				+ " complexes (before " + initialSize + ")");
		log.info("From the " + (initialSize - proteinComplexes.size()) + " discarded, " + big + " were too big (>"
				+ maxComplexSizeInReferenceSetForLearning + ") and " + small + " were too small (<"
				+ minComplexSizeInReferenceSetForLearning + ")");
	}

	public ProteinComplex getProteinComplexesByID(String complexID) {
		return proteinComplexesByID.get(complexID);
	}
}
