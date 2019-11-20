package edu.scripps.yates.pcomplex.distances;

import org.apache.commons.math3.ml.distance.EuclideanDistance;

import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.list.TDoubleList;

public class EuclideanDistanceCalculator {
	private static final EuclideanDistance euclideanDistance = new EuclideanDistance();

	public static double getEuclideanDistance(SeparationExperiment exp, String acc1, String acc2, DataType dataType) {
		final TDoubleList profile1 = exp.getDoubleNormalizedElutionProfile(acc1, dataType);
		final TDoubleList profile2 = exp.getDoubleNormalizedElutionProfile(acc2, dataType);
		final double euDistance = compute(profile1.toArray(), profile2.toArray());
		return euDistance;
	}

	public static double compute(double[] array1, double[] array2) {
		final double euDistance = euclideanDistance.compute(array1, array2);
		return euDistance;
	}
}
