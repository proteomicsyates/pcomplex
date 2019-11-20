package edu.scripps.yates.pcomplex.distances;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.list.array.TDoubleArrayList;

public class PearsonCorrelation {
	private static final PearsonsCorrelation pearson = new PearsonsCorrelation();

	public static double getPearsonCorrelationCoefficient(SeparationExperiment exp, String acc1, String acc2,
			DataType dataType) {
		final TDoubleArrayList profile1 = exp.getElutionProfile(acc1, dataType);
		final TDoubleArrayList profile2 = exp.getElutionProfile(acc2, dataType);

		final double correlation = correlationCoefficient(profile1.toArray(), profile2.toArray());

		return correlation;
	}

	public static double correlationCoefficient(double[] array1, double[] array2) {
		final double correlation = pearson.correlation(array1, array2);
		return correlation;
	}

	public static double getPearsonCorrelationPValue(SeparationExperiment exp, String acc1, String acc2,
			DataType dataType) {
		final double[][] rectangularData = getRectangularData(exp, acc1, acc2, dataType);
		final PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation(rectangularData);
		final RealMatrix matrix = pearsonsCorrelation.getCorrelationPValues();
		final double entry = matrix.getEntry(0, 1);
		return entry;
	}

	/**
	 * Get a matrix, which is two columns, one per each profile of each protein
	 * 
	 * @param exp
	 * @param acc1
	 * @param acc2
	 * @param dataType
	 * @return
	 */
	private static double[][] getRectangularData(SeparationExperiment exp, String acc1, String acc2,
			DataType dataType) {
		final TDoubleArrayList profile1 = exp.getElutionProfile(acc1, dataType);
		final TDoubleArrayList profile2 = exp.getElutionProfile(acc2, dataType);
		final double[][] ret = new double[profile1.size()][2];
		for (int i = 0; i < profile1.size(); i++) {
			final double d1 = profile1.get(i);
			final double d2 = profile2.get(i);
			ret[i][0] = d1;
			ret[i][1] = d2;
		}
		return ret;
	}
}
