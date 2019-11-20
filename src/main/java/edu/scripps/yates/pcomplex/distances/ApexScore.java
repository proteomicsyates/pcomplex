package edu.scripps.yates.pcomplex.distances;

import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.list.array.TDoubleArrayList;

public class ApexScore {

	/**
	 * Returns 1 if the two proteins have the maximum of its elution profiles at
	 * the same fraction.
	 * 
	 * @param exp
	 * @param protein1
	 * @param protein2
	 * @param dataType
	 * @return
	 */
	public static double getApexScore(SeparationExperiment exp, String protein1, String protein2, DataType dataType) {
		final TDoubleArrayList profile1 = exp.getElutionProfile(protein1, dataType);
		final TDoubleArrayList profile2 = exp.getElutionProfile(protein2, dataType);
		int index1 = -1;
		int index2 = -1;
		double max1 = -Double.MAX_VALUE;
		double max2 = -Double.MAX_VALUE;
		for (int index = 0; index < profile1.size(); index++) {
			if (max1 < profile1.get(index)) {
				max1 = profile1.get(index);
				index1 = index;
			}
			if (max2 < profile2.get(index)) {
				max2 = profile2.get(index);
				index2 = index;
			}
		}
		if (index1 == index2) {
			return 1.0;
		}
		return 0.0;
	}

}
