package edu.scripps.yates.pcomplex.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.sequence.PeptideSequenceProperties;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinComplex {
	private final static Logger log = Logger.getLogger(ProteinComplex.class);
	public static boolean useGeneToPrint = false;
	public static String separator = "-";

	private final Map<String, ProteinComponent> componentsByAccAndGene = new THashMap<String, ProteinComponent>();
	private String id;
	private String name;
	private final List<ProteinComponent> componentList = new ArrayList<ProteinComponent>();
	private final Set<ProteinComponent> componentSet = new THashSet<ProteinComponent>();
	private Double mw;
	private Double pi;
	private int hashcode = -1;
	private String key;
	private double density;
	private double significance;
	private final TObjectDoubleMap<String> maxOverlapsPerDB = new TObjectDoubleHashMap<String>();
	private boolean sorted;
	private boolean removed = false;
	private String organism;
	private List<String> componentsAccsAndGenes;
	private double abundance;
	private boolean known;
	private String goID;
	private String corumID;
	private String intActID;
	private boolean useIdInEqualsMethod;

	public ProteinComplex(String id) {
		this.id = id;
		parseID(id);
	}

	private void parseID(String id2) {
		final List<String> triplets = new ArrayList<String>();
		if (id2.contains(",")) {
			final String[] split = id2.split(",");
			for (int i = 0; i < split.length; i++) {
				if (split[i].contains(";")) {
					triplets.add(split[i]);
				} else {
					// in this case, it is because is a comma inside of a triplet, so add it to the
					// previous one
					String previousTriplet = triplets.get(triplets.size() - 1);
					previousTriplet += "," + split[i];
					triplets.set(triplets.size() - 1, previousTriplet);
				}

			}
		} else {
			triplets.add(id2);
		}
		for (final String triplet : triplets) {

			if (triplet.contains(";")) {
				final String[] split = triplet.split(";");
				final String idType = split[1];
				final String id = split[2];

				String name = null;
				if ("CORUM".equals(idType)) {
					corumID = id;
					name = split[3];
				} else if ("QuickGO".equals(idType)) {
					goID = id;
				} else if ("IntAct".equals(idType)) {
					intActID = id;
					name = split[3];
				}
			}
		}
	}

	public void setId(String id) {
		this.id = id;
	}

	public void addComponent(ProteinComponent component) {
		if (componentSet.contains(component)) {
			return;
		}
		if (component.getAcc() == null && "null".equals(component.getGene())) {
			log.info("asdf");
		}
		sorted = false;
		key = null;
		hashcode = -1;
		if (component.getAcc() != null) {
			componentsByAccAndGene.put(component.getAcc(), component);
		}
		if (component.getGene() != null) {
			componentsByAccAndGene.put(component.getGene(), component);
		}

		componentSet.add(component);
		// reset sorted
		// if (!componentList.contains(component)) {
		componentList.add(component);

		// }
	}

	public Map<String, ProteinComponent> getComponents() {
		return componentsByAccAndGene;
	}

	public List<ProteinComponent> getComponentList() {
		sortComponentList();
		return componentList;
	}

	public String getId() {
		return id;
	}

	/**
	 * Returns TRUE if all components of the complex are in the input set of
	 * proteinKeys (accs or gene names)
	 * 
	 * @param keys
	 * @return
	 */
	public boolean isFullyRepresented(Set<String> proteinKeys) {
		for (final String component : getComponentsAccAndGenes()) {
			if (!proteinKeys.contains(component)) {
				return false;
			}
		}
		return true;
	}

	private List<String> getComponentsAccAndGenes() {
		if (componentsAccsAndGenes == null) {
			componentsAccsAndGenes = new ArrayList<String>();
			componentsAccsAndGenes.addAll(componentsByAccAndGene.keySet());
		}
		return componentsAccsAndGenes;
	}

	/**
	 * Returns TRUE if at least ONE component of the complex is in the input set of
	 * proteinKeys (accs or gene names)
	 * 
	 * @param keys
	 * @return
	 */
	public boolean isPartiallyRepresented(Set<String> proteinKeys) {
		for (final String component : getComponentsAccAndGenes()) {
			if (proteinKeys.contains(component)) {
				return true;
			}
		}
		return false;
	}

	public void setName(String proteinComplexName) {
		name = proteinComplexName;
	}

	public String getName() {
		return name;
	}

	@Override
	public String toString() {
		return toString(useGeneToPrint, separator);
	}

	public String toString(boolean useGene, String separator) {
		final StringBuilder sb = new StringBuilder();
		if (getId() != null) {
			sb.append(getId() + " ");
		}
		if (getName() != null) {
			sb.append(getName());
		}
		boolean writeBraket = false;
		if (!"".equals(sb.toString())) {
			sb.append(": [");
			writeBraket = true;
		}

		int i = 0;
		for (final ProteinComponent component : getComponentList()) {
			if (i > 0) {
				sb.append(separator);
			}
			if (!useGene) {
				sb.append(component.getAcc());
			} else {
				sb.append(component.getGene());
			}
			i++;
		}
		if (writeBraket) {
			sb.append("]");
		}
		return sb.toString();
	}

	public List<ProteinComponent> getComponentsIn(Set<Protein> proteins) {
		final List<ProteinComponent> ret = new ArrayList<ProteinComponent>();
		final Set<String> inputs = new THashSet<String>();
		for (final Protein protein : proteins) {
			inputs.add(protein.getAcc());
			inputs.add(protein.getGene());
		}
		for (final ProteinComponent component : getComponentList()) {
			if (inputs.contains(component.getAcc()) || inputs.contains(component.getGene())) {
				ret.add(component);
			}
		}
		Collections.sort(ret, new Comparator<ProteinComponent>() {

			@Override
			public int compare(ProteinComponent o1, ProteinComponent o2) {
				return o1.getAcc().compareTo(o2.getAcc());
			}
		});
		return ret;

	}

	public List<ProteinComponent> getComponentsInProteinsAndGenes(Set<String> proteinsAndGenes) {
		final List<ProteinComponent> ret = new ArrayList<ProteinComponent>();

		for (final ProteinComponent component : getComponentList()) {

			if (proteinsAndGenes.contains(component.getAcc()) || proteinsAndGenes.contains(component.getGene())) {
				ret.add(component);
			} else if ((ProteinComplexAnalyzer.ORGANISMS != null && ProteinComplexAnalyzer.ORGANISMS.length > 0)
					|| getOrganism() != null) {
				// only if there is either acc or gene missing, otherwise it should be found
				if (component.getAcc() == null || component.getGene() == null) {
					UniprotGeneMapping geneMapping = null;
					if (ProteinComplexAnalyzer.ORGANISMS != null) {
						geneMapping = UniprotGeneMapping.getInstance(
								new File(ProteinComplexAnalyzer.uniprotReleasesFolder),
								ProteinComplexAnalyzer.ORGANISMS);
					} else {
						final String organism = getOrganism();
						geneMapping = UniprotGeneMapping
								.getInstance(new File(ProteinComplexAnalyzer.uniprotReleasesFolder), organism);
					}

					try {
						final Set<String> mapGeneToUniprotACC = geneMapping.mapGeneToUniprotACC(component.getGene());
						for (final String uniprot : mapGeneToUniprotACC) {
							if (proteinsAndGenes.contains(uniprot)) {
								ret.add(component);
							}
						}
					} catch (final IOException e) {
						e.printStackTrace();
					} catch (final IllegalArgumentException e2) {
					}
				}
			}
		}
		Collections.sort(ret, new Comparator<ProteinComponent>() {

			@Override
			public int compare(ProteinComponent o1, ProteinComponent o2) {
				if (o1.getAcc() != null && o2.getAcc() != null) {
					return o1.getAcc().compareTo(o2.getAcc());
				} else {
					return 1;
				}
			}
		});
		return ret;
	}

	public List<ProteinComponent> getComponentsNotIn(Set<Protein> proteins) {
		final List<ProteinComponent> ret = new ArrayList<ProteinComponent>();
		final Set<String> inputs = new THashSet<String>();
		for (final Protein protein : proteins) {
			inputs.add(protein.getAcc());
			inputs.add(protein.getGene());
		}
		for (final ProteinComponent component : getComponentList()) {
			if (!inputs.contains(component.getAcc()) && !inputs.contains(component.getGene())) {
				ret.add(component);
			}
		}
		Collections.sort(ret, new Comparator<ProteinComponent>() {

			@Override
			public int compare(ProteinComponent o1, ProteinComponent o2) {
				return o1.getAcc().compareTo(o2.getAcc());
			}
		});
		return ret;
	}

	public List<ProteinComponent> getComponentsNotInProteinsAndGenes(Set<String> proteinsAndGenes) {
		final List<ProteinComponent> ret = new ArrayList<ProteinComponent>();

		for (final ProteinComponent component : getComponentList()) {
			if (!proteinsAndGenes.contains(component.getAcc()) || !proteinsAndGenes.contains(component.getGene())) {
				ret.add(component);
			}
		}
		Collections.sort(ret, new Comparator<ProteinComponent>() {

			@Override
			public int compare(ProteinComponent o1, ProteinComponent o2) {
				final String acc1 = o1.getKey();
				final String acc2 = o2.getKey();
				return acc1.compareTo(acc2);
			}
		});
		return ret;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ProteinComplex) {
			final ProteinComplex complex = (ProteinComplex) obj;
			if (useIdInEqualsMethod) {
				return complex.getId().equals(getId()) && complex.getComponentsKey().equals(getComponentsKey());
			}
			return complex.getComponentsKey().equals(getComponentsKey());
		}
		return super.equals(obj);
	}

	public double getMW(UniprotGeneMapping geneMapping, UniprotProteinLocalRetriever uplr) throws IOException {
		if (mw == null) {
			final TDoubleArrayList allsubunitsMWs = new TDoubleArrayList();
			for (final ProteinComponent component : getComponentList()) {
				final Set<String> uniprotAccs = getUniprotACCsFromComponent(component, geneMapping, uplr);

				final TDoubleArrayList ambiguousEntriesMWs = new TDoubleArrayList();
				for (final String acc : uniprotAccs) {
					final Entry entry = uplr.getAnnotatedProtein(null, acc).get(acc);
					if (entry != null) {
						final Double mw = UniprotEntryUtil.getMolecularWeightInDalton(entry);
						if (mw != null) {
							ambiguousEntriesMWs.add(mw);
						}
					}
				}

				if (!ambiguousEntriesMWs.isEmpty()) {
					allsubunitsMWs.add(Maths.mean(ambiguousEntriesMWs));
				}
			}
			mw = allsubunitsMWs.sum();
		}
		return mw;
	}

	private Set<String> getUniprotACCsFromComponent(ProteinComponent component, UniprotGeneMapping geneMapping,
			UniprotProteinLocalRetriever uplr) throws IOException {

		Set<String> uniprotAccs = null;
		final String uniprotACC = FastaParser.getUniProtACC(component.getAcc());
		if (uniprotACC != null) {
			uniprotAccs = new THashSet<String>();
			uniprotAccs.add(uniprotACC);
		} else {
			final Set<String> mapGeneToUniprotACC = geneMapping.mapGeneToUniprotACC(component.getGene());
			uniprotAccs = new THashSet<String>();
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, mapGeneToUniprotACC);
			for (final String acc : mapGeneToUniprotACC) {
				if (annotatedProteins.containsKey(acc)) {
					final Entry entry = annotatedProteins.get(acc);

					if (UniprotEntryUtil.isSwissProt(entry)) {
						uniprotAccs.add(acc);
					}
				}
			}
			if (uniprotAccs.isEmpty()) {
				uniprotAccs.addAll(mapGeneToUniprotACC);
			}
		}
		return uniprotAccs;
	}

	public double getPI(UniprotGeneMapping geneMapping, UniprotProteinLocalRetriever uplr) throws IOException {
		if (pi == null) {
			final TDoubleArrayList allsubunitsPIs = new TDoubleArrayList();
			for (final ProteinComponent component : getComponentList()) {
				final Set<String> uniprotAccs = getUniprotACCsFromComponent(component, geneMapping, uplr);

				final TFloatArrayList ambiguousEntriesPIs = new TFloatArrayList();
				for (final String acc : uniprotAccs) {
					final Entry entry = uplr.getAnnotatedProtein(null, acc).get(acc);
					if (entry != null) {
						final Float pi = PeptideSequenceProperties
								.calculatepI(UniprotEntryUtil.getProteinSequence(entry));
						if (pi != null) {
							ambiguousEntriesPIs.add(pi);
						}
					}
				}

				if (!ambiguousEntriesPIs.isEmpty()) {
					allsubunitsPIs.add(Maths.mean(ambiguousEntriesPIs));
				}
			}
			pi = allsubunitsPIs.sum();
		}
		return pi;
	}

	public Set<ProteinComponent> getComponentSet() {
		return componentSet;
	}

	public List<String> getComponentListName(boolean useGeneToPrint) {
		final List<String> ret = new ArrayList<String>();
		for (final ProteinComponent component : getComponentList()) {
			final String componentName = component.toString(useGeneToPrint);
			ret.add(componentName);
		}
		return ret;
	}

	private void sortComponentList() {
		if (!sorted) {
			Collections.sort(componentList, new Comparator<ProteinComponent>() {

				@Override
				public int compare(ProteinComponent o1, ProteinComponent o2) {
					return o1.toString().compareTo(o2.toString());
				}
			});
			sorted = true;
		}
	}

	public List<String> getComponentListKey() {
		final List<String> ret = new ArrayList<String>();
		for (final ProteinComponent component : getComponentList()) {
			ret.add(component.getKey());
		}
		return ret;
	}

	public String getComponentsKey() {
		if (key == null) {
			final StringBuilder sb = new StringBuilder();
			for (final String key : getComponentListKey()) {
				if (!"".equals(sb.toString())) {
					sb.append("-");
				}
				sb.append(key);
			}
			key = sb.toString();
		}
		return key;
	}

	public String getComponentsGeneString(String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinComponent pc : getComponentList()) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			final String gene = pc.getGene();
			if (gene != null) {
				sb.append(gene);
			} else {
				sb.append(pc.getAcc());
			}
		}
		return sb.toString();
	}

	public String getComponentsProteinsAccString(String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinComponent pc : getComponentList()) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			final String acc = pc.getAcc();
			if (acc != null) {
				sb.append(acc);
			} else {
				sb.append(pc.getGene());
			}
		}
		return sb.toString();
	}

	@Override
	public int hashCode() {
		if (hashcode == -1) {
			String key = getComponentsKey();
			if (useIdInEqualsMethod) {
				key += getId();
			}
			hashcode = HashCodeBuilder.reflectionHashCode(key);
		}
		return hashcode;
	}

	public int size() {
		return componentList.size();
	}

	public void setDensity(double density) {
		this.density = density;
	}

	public double getDensity() {
		return density;
	}

	public void setSignificance(double significance) {
		this.significance = significance;
	}

	public double getSignificance() {
		return significance;
	}

	/**
	 * Gets the maximum overlap between this complex and a set of complexes in a
	 * {@link ProteinComplexDB}
	 * 
	 * @param db
	 * @return
	 */
	public double getMaxOverlap(ProteinComplexDB db) {
		if (!maxOverlapsPerDB.containsKey(db.getName())) {
			final double max = getMaxOverlap(db.getProteinComplexes());
			maxOverlapsPerDB.put(db.getName(), max);
		}
		return maxOverlapsPerDB.get(db.getName());
	}

	/**
	 * Gets the max overlap between this complex and a collection of complexes
	 * 
	 * @param complex2
	 * @return
	 */
	public double getMaxOverlap(Collection<ProteinComplex> complexes) {
		double max = -Double.MAX_VALUE;
		for (final ProteinComplex complex2 : complexes) {
			final double overlap = ClusterEvaluation.getOverlap(complex2, this);
			if (overlap > max) {
				max = overlap;
			}
		}
		return max;
	}

	/**
	 * Gets the overlap between this complex and complex2
	 * 
	 * @param complex2
	 * @return
	 */
	public double getOverlap(ProteinComplex complex2) {
		return ClusterEvaluation.getOverlap(this, complex2);
	}

	public boolean isRemoved() {
		return removed;
	}

	public void setRemoved(boolean removed) {
		this.removed = removed;
	}

	public String getOrganism() {
		return organism;
	}

	public void setOrganism(String organism) {
		this.organism = organism;
	}

	public void setAbundance(double abundance) {
		this.abundance = abundance;

	}

	public double getAbundance() {
		return abundance;
	}

	public void setKnown(boolean b) {
		known = b;

	}

	public boolean isKnown() {
		return known;
	}

	public boolean isSourceQuickGO() {
		return goID != null;
	}

	public boolean isSourceCorum() {
		return corumID != null;
	}

	public boolean isSourceIntAct() {
		return intActID != null;
	}

	public void setUseIdInEqualsMethod(boolean b) {
		useIdInEqualsMethod = b;
	}
}
