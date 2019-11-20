package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.DistanceMeasure;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import smile.netlib.NLMatrix;

public abstract class NaiveBayesClassification {
	private final static Logger log = Logger.getLogger(NaiveBayesClassification.class);
	protected final TrueClassifier trueClassifier;

	public NaiveBayesClassification(TrueClassifier trueClassifier) {
		this.trueClassifier = trueClassifier;
	}

	/**
	 * For each proteinPair,
	 * 
	 * @param matrixMap
	 * @param proteinList
	 * @return
	 * @throws IOException
	 */
	protected int[] getDataSetClasses(List<String> proteinList) throws IOException {
		final int numPairs = proteinList.size() * (proteinList.size() - 1) / 2;
		final int[] ret = new int[numPairs];
		int numPair = 0;
		for (int row = 0; row < proteinList.size(); row++) {
			final String protein1 = proteinList.get(row);
			for (int col = row + 1; col < proteinList.size(); col++) {
				final String protein2 = proteinList.get(col);
				final ClassLabel classLabel = trueClassifier.getClassLabel(protein1, protein2);
				ret[numPair] = classLabel.getClassNumber();
				numPair++;
			}
		}

		return ret;
	}

	/**
	 * For each protein pair, gets the values of all the matrices
	 * 
	 * @param matrixMap
	 * @param proteinList
	 * @return
	 */
	protected double[][] getDataSetData(Map<DistanceMeasure, NLMatrix> matrixMap) {
		log.info("Extracting dataset from matrix");
		final int ncols = matrixMap.values().iterator().next().ncols();
		final int nrows = matrixMap.values().iterator().next().nrows();
		final int numProteinPairs = ncols * (nrows - 1) / 2;
		final DistanceMeasure[] distances = DistanceMeasure.values();

		final double[][] ret = new double[numProteinPairs][distances.length];
		for (int dist = 0; dist < distances.length; dist++) {
			final DistanceMeasure distance = distances[dist];
			final NLMatrix matrix = matrixMap.get(distance);
			int numPair = 0;
			for (int row = 0; row < nrows; row++) {
				for (int col = row + 1; col < ncols; col++) {
					final double num = matrix.get(row, col);
					if (Double.isNaN(num)) {
						throw new IllegalArgumentException("This shoudn't happen");
						// TODO this shoudn't happen. Check why this number is
						// NaN
					}
					ret[numPair][dist] = num;
					numPair++;
				}
			}
		}
		return ret;
	}

	public abstract MyClassifier naiveBayesTraining() throws Exception;

}
