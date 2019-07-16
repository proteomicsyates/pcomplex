package edu.scripps.yates.pcomplex.cofractionation.training;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.DistanceMeasure;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.utilities.dates.DatesUtil;
import gnu.trove.map.hash.TIntIntHashMap;
import smile.classification.NaiveBayes;
import smile.netlib.NLMatrix;
import smile.stat.distribution.Distribution;
import smile.stat.distribution.GaussianMixture;
import smile.validation.CrossValidation;

public abstract class NaiveBayesSmile extends NaiveBayesClassification {
	private static final Logger log = Logger.getLogger(NaiveBayesSmile.class);

	public void naiveBayesTrainingOLD(Map<DistanceMeasure, NLMatrix> matrixMap, List<String> proteinList)
			throws IOException {
		log.info("Starting Naive Bayesian approach");
		final double[][] x = getDataSetData(matrixMap);
		final int[] y = getDataSetClasses(proteinList);
		//
		final int n = x.length;
		final int k = 15;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		final TIntIntHashMap labelFrequesncy = new TIntIntHashMap();
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final int numClasses = ClassLabel.values().length;
			final int inputDataDimension = DistanceMeasure.values().length;
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.MULTINOMIAL, numClasses, inputDataDimension);

			bayes.learn(trainx, trainy);
			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);
			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (labelFrequesncy.contains(label)) {
					labelFrequesncy.put(label, labelFrequesncy.get(label) + 1);
				} else {
					labelFrequesncy.put(label, 1);
				}
				if (label != -1) {
					// System.out.println("Pair: " + testx[j][0] + "," +
					// testx[j][1] + "," + testx[j][2] + ","
					// + testx[j][3] + "," + testx[j][4] + "\t" +
					// ClassLabel.values()[label].name()
					// + " and correct is " +
					// ClassLabel.values()[testy[j]].name());
					total++;
					if (testy[j] != label) {
						error++;
					}
					final double[] posteriori = new double[ClassLabel.values().length];
					bayes.predict(testx[j], posteriori);
					int q = 0;
					for (final double d : posteriori) {
						System.out.println(ClassLabel.values()[q++] + ": " + d);
					}
				} else {
					// log.info("why: " + label);

				}
			}
			final double[] priori = bayes.getPriori();
			for (final double d : priori) {
				System.out.println(d);
			}
		}
		for (final int i : labelFrequesncy.keys()) {
			System.out.println("Num pairs with class: " + i + ": " + labelFrequesncy.get(i));
		}
		System.out.format("Batch Multinomial error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);
	}

	@Override
	public ClassificationResult naiveBayesTraining(Map<DistanceMeasure, NLMatrix> matrixMap, List<String> proteinList)
			throws IOException {
		log.info("Starting Naive Bayesian approach");
		final double[][] x = getDataSetData(matrixMap);
		final int[] y = getDataSetClasses(proteinList);
		//
		final int n = x.length;
		final int kk = 15;
		final CrossValidation cv = new CrossValidation(n, kk);
		int error = 0;
		int total = 0;
		final TIntIntHashMap labelFrequesncy = new TIntIntHashMap();
		for (int l = 0; l < kk; l++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[l]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[l]);

			final int p = trainx[0].length;
			final int k = smile.math.Math.max(trainy) + 1;

			final double[] priori = new double[k];
			final Distribution[][] condprob = new Distribution[k][p];
			for (int i = 0; i < k; i++) {
				priori[i] = 1.0 / k;
				for (int j = 0; j < p; j++) {
					final ArrayList<Double> axi = new ArrayList<>();
					for (int m = 0; m < trainx.length; m++) {
						if (trainy[m] == i) {
							axi.add(trainx[m][j]);
						}
					}

					final double[] xi = new double[axi.size()];
					for (int m = 0; m < xi.length; m++) {
						xi[m] = axi.get(m);
					}
					log.info("Calculating conditional probability of class " + i + "(" + ClassLabel.values()[i].name()
							+ ") and measurement " + j + " is done. Using " + xi.length + " data points");
					final long t0 = System.currentTimeMillis();
					condprob[i][j] = new GaussianMixture(xi, k);
					final long t1 = System.currentTimeMillis();
					log.info("Conditional probability of class " + i + "(" + ClassLabel.values()[i].name()
							+ ") and measurement " + j + " is done in "
							+ DatesUtil.getDescriptiveTimeFromMillisecs((t1 - t0)));
				}
			}
			final NaiveBayes bayes = new NaiveBayes(priori, condprob);
			final double[] posteriori = new double[ClassLabel.values().length];
			final double[][] testx = smile.math.Math.slice(x, cv.test[l]);
			final int[] testy = smile.math.Math.slice(y, cv.test[l]);
			for (int j = 0; j < testx.length; j++) {

				final int label = bayes.predict(x[testy[j]], posteriori);
				if (labelFrequesncy.contains(label)) {
					labelFrequesncy.put(label, labelFrequesncy.get(label) + 1);
				} else {
					labelFrequesncy.put(label, 1);
				}
				if (label != -1) {
					// System.out.println("Pair: " + testx[j][0] + "," +
					// testx[j][1] + "," + testx[j][2] + ","
					// + testx[j][3] + "," + testx[j][4] + "\t" +
					// ClassLabel.values()[label].name()
					// + " and correct is " +
					// ClassLabel.values()[testy[j]].name());
					total++;
					if (y[cv.test[l][j]] != label) {
						error++;
					}

					int q = 0;
					for (final double d : posteriori) {
						System.out.println(ClassLabel.values()[q++].name() + ": " + d);
					}
				} else {
					// log.info("why: " + label);

				}
			}
			for (final double d : priori) {
				System.out.println(d);
			}
		}
		for (final int i : labelFrequesncy.keys()) {
			System.out.println("Num pairs with class: " + i + ": " + labelFrequesncy.get(i));
		}
		System.out.format("Batch Multinomial error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);

		// TODO
		return null;
	}

	@Override
	public abstract ClassLabel getClassLabel(String protein, String protein2) throws IOException;
}
