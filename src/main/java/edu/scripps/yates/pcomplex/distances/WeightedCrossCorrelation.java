package edu.scripps.yates.pcomplex.distances;

import org.geneontology.oboedit.datamodel.Datatype;

import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

public class WeightedCrossCorrelation {

	/**
	 * Calculates the weighted cross correlation from the raw elution profiles
	 * (of {@link Datatype}) of 2 proteins in a {@link SeparationExperiment}
	 * 
	 * @param experiment
	 * @param protein1
	 * @param protein2
	 * @param dataType
	 * @return
	 */
	public static double getWeightedCrossCorrelation(SeparationExperiment experiment, String protein1, String protein2,
			DataType dataType) {
		final TDoubleArrayList profile1 = experiment.getElutionProfile(protein1, dataType);
		final TDoubleArrayList profile2 = experiment.getElutionProfile(protein2, dataType);
		final double wcc12 = weigthedCrossCorrelation(profile1, profile2);
		final double wcc11 = weigthedCrossCorrelation(profile1, profile1);
		final double wcc22 = weigthedCrossCorrelation(profile2, profile2);
		if (wcc11 == 0.0 && wcc22 == 0.0) {
			return 0.0;
		}
		final double ret = wcc12 / Math.sqrt(wcc11 * wcc22);
		return ret;
	}

	public static double weigthedCrossCorrelation(TDoubleList profile1, TDoubleList profile2) {
		double ret = 0.0;
		for (int offset = -1; offset <= 1; offset++) {
			final TDoubleList subProfile1 = profile1.subList(Math.max(0, offset),
					Math.min(profile1.size(), profile1.size() + offset));
			final TDoubleList subProfile2 = profile2.subList(Math.max(0, -offset),
					Math.min(profile2.size(), profile2.size() - offset));
			final double cc = crossCorrelation(subProfile1.toArray(), subProfile2.toArray());
			final double weight = triangleFunction(offset);
			ret += weight * cc;
		}
		return ret;
	}

	public static double triangleFunction(double x) {
		if (Math.abs(x) >= 2.0) {
			return 0;
		}
		return 1 - Math.abs(x) / 2.0;
	}

	public static double crossCorrelation(double[] a, double b[]) {
		if (a.length != b.length) {
			throw new IllegalArgumentException("Arrays must have same length");
		}
		double ret = 0.0;
		for (int index = 0; index < a.length; index++) {
			ret += a[index] * b[index];
		}
		return ret;
	}

}
