package edu.scripps.yates.pcomplex.model;

import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class Fraction {
	private static final Logger log = Logger.getLogger(Fraction.class);
	private final String name;
	private final int number;
	private final Set<Protein> proteins = new THashSet<Protein>();
	private final Map<String, Protein> proteinsByAcc = new THashMap<String, Protein>();
	private final Map<String, Set<Protein>> proteinsByGene = new THashMap<String, Set<Protein>>();
	private THashSet<String> proteinsAndGenes;

	public Fraction(String name, int number) {
		this.name = name;
		this.number = number;
	}

	public void addProtein(Protein protein) {
		proteins.add(protein);
		if (!proteinsByAcc.containsKey(protein.getAcc())) {
			proteinsByAcc.put(protein.getAcc(), protein);
		} else {
			log.warn("Protein '" + protein.getAcc() + "' already present in fraction " + this);
		}
		if (proteinsByGene.containsKey(protein.getGene())) {
			proteinsByGene.get(protein.getGene()).add(protein);
		} else {
			final Set<Protein> set = new THashSet<Protein>();
			set.add(protein);
			proteinsByGene.put(protein.getGene(), set);
		}
	}

	@Override
	public String toString() {
		return number + "-" + name;
	}

	public String getName() {
		return name;
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

}
