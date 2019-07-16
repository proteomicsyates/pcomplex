package edu.scripps.yates.pcomplex.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.builder.HashCodeBuilder;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.sequence.PeptideSequenceProperties;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.set.hash.THashSet;

public class ProteinComplex {
	private final Set<String> components = new THashSet<String>();
	private final String id;
	private String name;
	private ArrayList<String> sortedComponents;
	private final boolean geneNames;
	private ArrayList<String> componentList;
	private Double mw;
	private Double pi;

	private ProteinComplex(String id) {
		this(id, false);
	}

	public ProteinComplex(String id, boolean geneNames) {
		this.id = id;
		this.geneNames = geneNames;
	}

	public void addComponent(String component) {
		components.add(component);
		// reset sorted
		sortedComponents = null;
	}

	public Set<String> getComponents() {
		return components;
	}

	public List<String> getComponentList() {
		if (componentList == null) {
			componentList = new ArrayList<String>();
			componentList.addAll(getComponents());
			Collections.sort(componentList);
		}
		return componentList;
	}

	public boolean contains(Protein protein) {
		if (components.contains(protein.getAcc()) || components.contains(protein.getGene())) {
			return true;
		}
		return false;
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
		for (final String component : components) {
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
		for (final String component : components) {
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
		final StringBuilder sb = new StringBuilder();
		sb.append(getId() + "\t");
		if (getName() != null) {
			sb.append(getName() + "\n");
		} else {
			sb.append("no-name\n");
		}
		for (final String component : getSortedComponents()) {
			sb.append(component + "\t");
		}
		return sb.toString();
	}

	private List<String> getSortedComponents() {
		if (sortedComponents == null) {
			sortedComponents = new ArrayList<String>();
			sortedComponents.addAll(components);
			Collections.sort(sortedComponents);
		}
		return sortedComponents;
	}

	public List<String> getComponentsIn(Set<Protein> proteins) {
		final List<String> ret = new ArrayList<String>();
		final Set<String> inputs = new THashSet<String>();
		for (final Protein protein : proteins) {
			inputs.add(protein.getAcc());
			inputs.add(protein.getGene());
		}
		for (final String component : components) {
			if (inputs.contains(component)) {
				ret.add(component);
			}
		}
		Collections.sort(ret);
		return ret;
	}

	public List<String> getComponentsInProteinsAndGenes(Set<String> proteinsAndGenes) {
		final List<String> ret = new ArrayList<String>();

		for (final String component : components) {

			if (proteinsAndGenes.contains(component)) {
				ret.add(component);
			} else {
				final UniprotGeneMapping geneMapping = UniprotGeneMapping.getInstance(
						new File(ProteinComplexAnalyzer.uniprotReleasesFolder), ProteinComplexAnalyzer.TAXONOMY);
				try {
					final Set<String> mapGeneToUniprotACC = geneMapping.mapGeneToUniprotACC(component);
					for (final String uniprot : mapGeneToUniprotACC) {
						if (proteinsAndGenes.contains(uniprot)) {
							ret.add(component);
						}
					}
				} catch (final IOException e) {
					e.printStackTrace();
				}
			}
		}
		Collections.sort(ret);
		return ret;
	}

	public List<String> getComponentsNotIn(Set<Protein> proteins) {
		final List<String> ret = new ArrayList<String>();
		final Set<String> inputs = new THashSet<String>();
		for (final Protein protein : proteins) {
			inputs.add(protein.getAcc());
			inputs.add(protein.getGene());
		}
		for (final String component : components) {
			if (!inputs.contains(component)) {
				ret.add(component);
			}
		}
		Collections.sort(ret);
		return ret;
	}

	public List<String> getComponentsNotInProteinsAndGenes(Set<String> proteinsAndGenes) {
		final List<String> ret = new ArrayList<String>();

		for (final String component : components) {
			if (!proteinsAndGenes.contains(component)) {
				ret.add(component);
			}
		}
		Collections.sort(ret);
		return ret;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ProteinComplex) {
			final ProteinComplex complex = (ProteinComplex) obj;
			if (complex.getComponents().size() != getComponents().size()) {
				return false;
			}
			for (final String string : components) {
				if (!complex.getComponents().contains(string)) {
					return false;
				}
			}
			for (final String string : complex.getComponents()) {
				if (!getComponents().contains(string)) {
					return false;
				}
			}
			return true;
		}
		return super.equals(obj);
	}

	@Override
	public int hashCode() {
		return HashCodeBuilder.reflectionHashCode(getSortedComponents());
	}

	public boolean isGeneNames() {
		return geneNames;
	}

	public double getMW(UniprotGeneMapping geneMapping, UniprotProteinLocalRetriever uplr) throws IOException {
		if (mw == null) {
			final TDoubleArrayList allsubunitsMWs = new TDoubleArrayList();
			for (final String componentKey : getComponentList()) {
				final Set<String> uniprotAccs = getUniprotACCsFromComponentKey(componentKey, geneMapping, uplr);

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

	private Set<String> getUniprotACCsFromComponentKey(String componentKey, UniprotGeneMapping geneMapping,
			UniprotProteinLocalRetriever uplr) throws IOException {
		Set<String> uniprotAccs = null;
		final String uniprotACC = FastaParser.getUniProtACC(componentKey);
		if (uniprotACC != null) {
			uniprotAccs = new THashSet<String>();
			uniprotAccs.add(uniprotACC);
		} else {
			final Set<String> mapGeneToUniprotACC = geneMapping.mapGeneToUniprotACC(componentKey);
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
			for (final String componentKey : getComponentList()) {
				final Set<String> uniprotAccs = getUniprotACCsFromComponentKey(componentKey, geneMapping, uplr);

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
}
