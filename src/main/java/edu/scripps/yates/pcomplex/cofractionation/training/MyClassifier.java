package edu.scripps.yates.pcomplex.cofractionation.training;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import weka.core.Instance;

public abstract class MyClassifier {

	public abstract ClassLabel classifyInstance(Instance instance) throws Exception;

	public abstract double[] distributionForInstance(Instance instance) throws Exception;

	public ClassificationResult evaluateClassifier(Dataset dataset) throws Exception {
		final ClassificationResult result = new ClassificationResult(this);
		int instanceIndex = 0;
		for (final String proteinPair : dataset.getProteinPairs()) {

			final ClassLabel trueClassLabel = dataset.getTrueClass(proteinPair);
			final Instance instance = dataset.getInstances().get(instanceIndex);
			final double[] probabilities = distributionForInstance(instance);
			final ClassLabel classifyInstance = classifyInstance(instance);
			// if (classifyInstance == 0.0) {
			// log.info("ASDF");
			// }

			final String protein1 = dataset.getFirstProtein(proteinPair);
			final String protein2 = dataset.getSecondProtein(proteinPair);
			final ProteinPairInteraction ppi = new ProteinPairInteraction(protein1, protein2, probabilities,
					trueClassLabel);
			result.addProteinInteraction(ppi);
			instanceIndex++;

		}
		return result;
	}

}
