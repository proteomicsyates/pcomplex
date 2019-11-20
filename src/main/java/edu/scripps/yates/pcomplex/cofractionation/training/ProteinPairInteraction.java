package edu.scripps.yates.pcomplex.cofractionation.training;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;

public class ProteinPairInteraction {
	private final String protein1Acc;
	private final String protein2Acc;
	private final double[] probabilities;
	private final ClassLabel trueLabel;
	private double precision;
	private ClassLabel predictedClass;
	private static final ClassLabel[] classLabels = ClassLabel.valuesArray();

	public ProteinPairInteraction(String protein1Acc, String protein2Acc, double[] probabilities,
			ClassLabel trueLabel) {
		this.protein1Acc = protein1Acc;
		this.protein2Acc = protein2Acc;
		this.probabilities = probabilities;
		this.trueLabel = trueLabel;
	}

	public String getProtein1Acc() {
		return protein1Acc;
	}

	public String getProtein2Acc() {
		return protein2Acc;
	}

	public double[] getProbabilities() {
		return probabilities;
	}

	public ClassLabel getTrueLabel() {
		return trueLabel;
	}

	public ClassLabel getPredictedLabel() {
		if (predictedClass == null) {
			double max = 0.0;
			int maxIndex = 0;
			for (int i = 0; i < probabilities.length; i++) {
				final double prob = probabilities[i];
				if (prob > max) {
					maxIndex = i;
					max = prob;
				}
			}
			if (max > 0) {

				final ClassLabel ret = classLabels[maxIndex];
				predictedClass = ret;
			}
		}
		return predictedClass;
	}

	public double getProbability(ClassLabel classLabel) {
		if (classLabel != null && classLabel.ordinal() < probabilities.length) {
			return probabilities[classLabel.ordinal()];
		}
		return Double.NaN;
	}

	public double getPrecision() {
		return precision;
	}

	public void setPrecision(double precision) {
		this.precision = precision;
	}
}
