package edu.scripps.yates.pcomplex.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class Fraction {
	private static final Logger log = Logger.getLogger(Fraction.class);

	private final int number;
	private final Set<Protein> proteins = new THashSet<Protein>();
	private final Map<String, Protein> proteinsByAcc = new THashMap<String, Protein>();
	private final Map<String, Set<Protein>> proteinsByGene = new THashMap<String, Set<Protein>>();
	private THashSet<String> proteinsAndGenes;
	private final THashSet<String> fractionNames = new THashSet<String>();

	private boolean containsReplicates = false;

	private final TIntList replicateNumbers = new TIntArrayList();

	public Fraction(String name, int number) {
		replicateNumbers.add(1);
		addFractionName(name);
		this.number = number;
	}

	public void addProtein(Protein protein) {
		proteins.add(protein);
		if (!proteinsByAcc.containsKey(protein.getAcc())) {
			proteinsByAcc.put(protein.getAcc(), protein);

			if (proteinsByGene.containsKey(protein.getGene())) {
				proteinsByGene.get(protein.getGene()).add(protein);
			} else {
				final Set<Protein> set = new THashSet<Protein>();
				set.add(protein);
				proteinsByGene.put(protein.getGene(), set);
			}
		} else {
			log.debug("Protein '" + protein.getAcc() + "' already present in fraction " + this);
			final int repNumber = proteinsByAcc.get(protein.getAcc()).addSpc(protein.getSpc());
			containsReplicates = true;
			if (!replicateNumbers.contains(repNumber)) {
				replicateNumbers.add(repNumber);
			}
		}

	}

	@Override
	public String toString() {
		return number + "-" + getFractionNamesString();
	}

	public String getFractionNamesString() {
		final List<String> ret = new ArrayList<String>();
		ret.addAll(fractionNames);
		Collections.sort(ret);
		final StringBuilder sb = new StringBuilder();
		for (final String fractionName : ret) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(fractionName);
		}
		return sb.toString();
	}

	public Set<Protein> getProteins() {
		return proteins;
	}

	public Set<String> getProteinsAndGenes() {
		if (proteinsAndGenes == null) {
			proteinsAndGenes = new THashSet<String>();
			for (final Protein protein : proteins) {
				proteinsAndGenes.add(protein.getAcc());
				if (protein.getGene() != null) {
					proteinsAndGenes.add(protein.getGene());
				}
			}
		}
		return proteinsAndGenes;
	}

	public Protein getProteinByAcc(String acc) {
		return proteinsByAcc.get(acc);
	}

	public Set<Protein> getProteinsByGene(String gene) {
		return proteinsByGene.get(gene);
	}

	public int getFractionNumber() {
		return number;
	}

	public Set<Protein> getProteinByComponent(ProteinComponent component) {
		final Set<Protein> ret = new THashSet<Protein>();
		if (component.getAcc() != null) {
			if (proteinsByAcc.containsKey(component.getAcc())) {
				ret.add(proteinsByAcc.get(component.getAcc()));
			}
		}
		if (component.getGene() != null) {
			if (proteinsByGene.containsKey(component.getGene())) {
				ret.addAll(proteinsByGene.get(component.getGene()));
			}
		}
		return ret;
	}

	public Set<Protein> getProteinByKey(String componentKey) {
		final Set<Protein> ret = new THashSet<Protein>();
		if (proteinsByAcc.containsKey(componentKey)) {
			ret.add(proteinsByAcc.get(componentKey));
		}
		if (proteinsByGene.containsKey(componentKey)) {
			ret.addAll(proteinsByGene.get(componentKey));
		}
		return ret;
	}

	public Set<ProteinComplex> getCompleteComplexes(ProteinComplexDB proteinComplexDB, int minNumComponentsInComplex) {
		return proteinComplexDB.getCompleteProteinComplexes(getProteinsAndGenes(), minNumComponentsInComplex);

	}

	public Set<ProteinComplex> getCompleteComplexes(ProteinComplexDB proteinComplexDB,
			ProteinComplexExistenceCriteria existenceCriteria, int minNumComponentsInComplex) {
		return proteinComplexDB.getCompleteProteinComplexes(getProteinsAndGenes(), existenceCriteria,
				minNumComponentsInComplex);

	}

	public void addFractionName(String fractionName) {
		this.fractionNames.add(fractionName);
	}

	public Set<String> getFractionNames() {
		return fractionNames;
	}

	public boolean containsReplicates() {

		return containsReplicates;
	}

	public TIntList getReplicates() {
		return replicateNumbers;
	}

}
