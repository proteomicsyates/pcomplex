package edu.scripps.yates.pcomplex.cofractionation.training;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import weka.classifiers.Classifier;
import weka.core.Instance;

public class MyClassifierWeka extends MyClassifier {
	private final Classifier wekaClassifier;
	private static final ClassLabel[] classes = ClassLabel.values();

	public MyClassifierWeka(Classifier classifier) {
		wekaClassifier = classifier;
	}

	@Override
	public ClassLabel classifyInstance(Instance instance) throws Exception {
		final double[] probabilities = wekaClassifier.distributionForInstance(instance);
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
			final ClassLabel ret = classes[maxIndex];
			return ret;
		}
		return null;

	}

	@Override
	public double[] distributionForInstance(Instance instance) throws Exception {
		return wekaClassifier.distributionForInstance(instance);
	}
}
