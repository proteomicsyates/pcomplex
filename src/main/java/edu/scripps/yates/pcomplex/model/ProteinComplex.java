package edu.scripps.yates.pcomplex.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.builder.HashCodeBuilder;

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
	public static boolean useGeneToPrint = false;
	public static String separator = "-";

	private final Map<String, ProteinComponent> componentsByAccAndGene = new THashMap<String, ProteinComponent>();
	private final String id;
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

	public ProteinComplex(String id) {
		this.id = id;

	}

	public void addComponent(ProteinComponent component) {
		if (componentSet.contains(component)) {
			return;
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
		for (final String component : componentsByAccAndGene.keySet()) {
			if (!proteinKeys.contains(component)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns TRUE if at least ONE component of the complex is in the input set of
	 * proteinKeys (accs or gene names)
	 * 
	 * @param keys
	 * @return
	 */
	public boolean isPartiallyRepresented(Set<String> proteinKeys) {
		for (final String component : componentsByAccAndGene.keySet()) {
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
			} else if (ProteinComplexAnalyzer.ORGANISM != null || getOrganism() != null) {
				final String organism = ProteinComplexAnalyzer.ORGANISM != null ? ProteinComplexAnalyzer.ORGANISM
						: getOrganism();
				final UniprotGeneMapping geneMapping = UniprotGeneMapping
						.getInstance(new File(ProteinComplexAnalyzer.uniprotReleasesFolder), organism);
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
				return o1.getAcc().compareTo(o2.getAcc());
			}
		});
		return ret;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ProteinComplex) {
			final ProteinComplex complex = (ProteinComplex) obj;

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
			ret.add(component.toString(useGeneToPrint));
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
			sb.append(pc.getGene());
		}
		return sb.toString();
	}

	public String getComponentsProteinsAccString(String separator) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinComponent pc : getComponentList()) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(pc.getAcc());
		}
		return sb.toString();
	}

	@Override
	public int hashCode() {
		if (hashcode == -1) {
			hashcode = HashCodeBuilder.reflectionHashCode(getComponentsKey());
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

	public double getMaxOverlap(ProteinComplexDB db) {
		if (!maxOverlapsPerDB.containsKey(db.getName())) {
			double max = -Double.MAX_VALUE;
			for (final ProteinComplex complex : db.getProteinComplexes()) {
				final double overlap = ClusterEvaluation.getOverlap(complex, this);
				if (overlap > max) {
					max = overlap;
				}
			}
			maxOverlapsPerDB.put(db.getName(), max);
		}
		return maxOverlapsPerDB.get(db.getName());
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

}
