package edu.scripps.yates.pcomplex.cofractionation.training;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;

public class ClassificationResult {

	private final double[][][] matrix;
	private final ClassLabel[][] trueLabelsMatrix;

	public ClassificationResult(int numRows, int numCols, int numClasses) {
		matrix = new double[numRows][numCols][numClasses];
		trueLabelsMatrix = new ClassLabel[numRows][numCols];
	}

	public void set(int i, int j, double[] probabilities, ClassLabel trueLabel) {
		for (int k = 0; k < probabilities.length; k++) {
			matrix[i][j][k] = probabilities[k];
		}
		trueLabelsMatrix[i][j] = trueLabel;
	}

}
