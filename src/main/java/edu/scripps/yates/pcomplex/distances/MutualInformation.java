package edu.scripps.yates.pcomplex.distances;

import java.util.List;
import java.util.Set;

import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class MutualInformation {
	/**
	 * Gets the mutual information between 2 proteins in a
	 * {@link SeparationExperiment}
	 * 
	 * @param exp
	 * @param protein1
	 * @param protein2
	 * @return
	 */
	public static double getMutualInformation(SeparationExperiment exp, String protein1, String protein2) {

		final double h12 = h(exp, protein1, protein2);
		final double h1 = h(exp, protein1);
		final double h2 = h(exp, protein2);
		final double mi = h12 - h1 - h2;
		return mi;
	}

	/**
	 * Entropy of protein1
	 * 
	 * @param fractions1
	 * @return
	 */
	private static double h(SeparationExperiment exp, String protein) {

		double ret = 0.0;
		for (int i = 0; i <= 1; i++) {
			final double p = p(exp, protein, i);
			if (p == 0.0) {
				continue;
			}
			ret += p * Maths.log(p, 2);
		}
		ret = -ret;
		return ret;
	}

	/**
	 * The ratio of the number of fractions in which a protein is present divided by
	 * the total number of fractions in the experiment
	 * 
	 * @param exp
	 * @param protein
	 * @param present if 1, means that we count the number of fractions in which a
	 *                protein is present. If 0, we count the number of fractions in
	 *                which a protein is NOT present.
	 */
	private static double p(SeparationExperiment exp, String protein, int present) {
		final Set<Integer> fractions = exp.getFractionsInWhichProteinIsPresent(protein);
		final int numPresent = fractions.size();
		final int totalFractions = exp.getFractions().size();
		if (present == 1) {
			final double p = 1.0 * numPresent / totalFractions;
			return p;
		} else {
			final int numNonPresent = totalFractions - numPresent;
			final double p = 1.0 * numNonPresent / totalFractions;
			return p;
		}
	}

	/**
	 * The ratio of the number of fractions in which a protein is present divided by
	 * the total number of fractions in the experiment
	 * 
	 * @param exp
	 * @param protein1
	 * @param protein2
	 * @param present  if 1, means that we count the number of fractions in which a
	 *                 protein is present. If 0, we count the number of fractions in
	 *                 which a protein is NOT present.
	 */
	private static double p(SeparationExperiment exp, String protein1, int present1, String protein2, int present2) {
		final Set<Integer> fractions1 = exp.getFractionsInWhichProteinIsPresent(protein1);
		final Set<Integer> fractions2 = exp.getFractionsInWhichProteinIsPresent(protein2);
		final TIntObjectHashMap<Set<Integer>> binaryFractions1 = getBinaryFractions(fractions1, exp.getFractions());
		final TIntObjectHashMap<Set<Integer>> binaryFractions2 = getBinaryFractions(fractions2, exp.getFractions());
		final Set<Integer> fractionsToConsider1 = binaryFractions1.get(present1);
		final Set<Integer> fractionsToConsider2 = binaryFractions2.get(present2);

		final Set<Integer> intersection = intersection(fractionsToConsider1, fractionsToConsider2);
		final int totalFractions = exp.getFractions().size();
		final double p = 1.0 * intersection.size() / totalFractions;
		return p;

	}

	/**
	 * Gets a map in which the key is 1 or 0, for the fractions in which a protein
	 * is present or is not present respectively
	 * 
	 * @param proteinFractions fractions in which a protein is present
	 * @param totalFractions
	 * @return
	 */
	private static TIntObjectHashMap<Set<Integer>> getBinaryFractions(Set<Integer> proteinFractions,
			List<Fraction> totalFractions) {
		final TIntObjectHashMap<Set<Integer>> ret = new TIntObjectHashMap<Set<Integer>>();
		ret.put(0, new THashSet<Integer>());
		ret.put(1, new THashSet<Integer>());
		for (final Fraction fraction : totalFractions) {
			if (proteinFractions.contains(fraction.getFractionNumber())) {
				ret.get(1).add(fraction.getFractionNumber());
			} else {
				ret.get(0).add(fraction.getFractionNumber());
			}
		}
		return ret;
	}

	private static Set<Integer> intersection(Set<Integer> fractionsToConsider1, Set<Integer> fractionsToConsider2) {
		final Set<Integer> intersection = new THashSet<Integer>();
		intersection.addAll(fractionsToConsider1);
		intersection.retainAll(fractionsToConsider2);
		return intersection;
	}

	/**
	 * joint entropy
	 * 
	 * @param fractions1
	 * @param fractions2
	 * @return
	 */
	private static double h(SeparationExperiment exp, String protein1, String protein2) {
		double ret = 0.0;
		for (int i = 0; i <= 1; i++) {
			for (int j = 0; j <= 1; j++) {
				final double p = p(exp, protein1, i, protein2, j);
				if (p == 0.0) {
					continue;
				}
				ret += p * Maths.log(p, 2);
			}
		}
		ret = -ret;
		return ret;
	}
}
