package edu.scripps.yates.pcomplex.distances;

import java.util.Set;

import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import smile.math.distance.JaccardDistance;

public class JaccardScore {
	/**
	 * Number of fractions that contain both proteins divided by the number of
	 * fractions that have at least one of the two proteins.
	 * 
	 * @param experiment2
	 * @param protein1
	 * @param protein2
	 * @param dataType
	 * @return
	 */
	public static double getJaccardScore(SeparationExperiment exp, String protein1, String protein2) {
		final Set<Integer> fractions1 = exp.getFractionsInWhichProteinIsPresent(protein1);
		final Set<Integer> fractions2 = exp.getFractionsInWhichProteinIsPresent(protein2);

		final double distance = 1 - JaccardDistance.d(fractions1, fractions2);

		return distance;
	}

}
