package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.IOException;

import org.apache.log4j.Logger;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.NaiveBayes;
import weka.core.Instances;

public class NaiveBayesWeka extends NaiveBayesClassification {
	private final static Logger log = Logger.getLogger(NaiveBayesWeka.class);
	private final int kFold;
	private final Instances dataSetTraining;

	public NaiveBayesWeka(int kFoldCrossValidation, double minCorrelation, Dataset dataset) throws IOException {
		super(dataset.getTrueClassifier());
		kFold = kFoldCrossValidation;

		dataSetTraining = dataset.getInstances();
	}

	@Override
	public MyClassifier naiveBayesTraining() throws Exception {
		log.info("Starting Naive Bayesian approach");

		final Classifier classifier = new NaiveBayes();

		for (int foldNumber = 0; foldNumber < kFold; foldNumber++) {
			final Instances trainCV = dataSetTraining.trainCV(kFold, foldNumber);

			classifier.buildClassifier(trainCV);

			final Instances testCV = dataSetTraining.testCV(kFold, foldNumber);
			final Evaluation eval = new Evaluation(trainCV);
			eval.evaluateModel(classifier, testCV);

			System.out.println("** Naive Bayes Evaluation with Datasets **");
			System.out.println(eval.toSummaryString());
			System.out.print(" the expression for the input data as per alogorithm is ");
			System.out.println(classifier);

		}

		return new MyClassifierWeka(classifier);

	}

}
