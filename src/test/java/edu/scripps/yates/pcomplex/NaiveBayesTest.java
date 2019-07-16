package edu.scripps.yates.pcomplex;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import smile.classification.NaiveBayes;
import smile.data.AttributeDataset;
import smile.data.parser.ArffParser;
import smile.feature.Bag;
import smile.stat.distribution.Distribution;
import smile.stat.distribution.GaussianMixture;
import smile.validation.CrossValidation;
import smile.validation.LOOCV;

public class NaiveBayesTest {

	String[] feature = { "outstanding", "wonderfully", "wasted", "lame", "awful", "poorly", "ridiculous", "waste",
			"worst", "bland", "unfunny", "stupid", "dull", "fantastic", "laughable", "mess", "pointless", "terrific",
			"memorable", "superb", "boring", "badly", "subtle", "terrible", "excellent", "perfectly", "masterpiece",
			"realistic", "flaws" };
	double[][] moviex;
	int[] moviey;

	public NaiveBayesTest() {
		final String[][] x = new String[2000][];
		final int[] y = new int[2000];

		try (final BufferedReader input = new BufferedReader(new InputStreamReader(
				new FileInputStream(new ClassPathResource("movie.txt").getFile().getAbsolutePath())))) {
			for (int i = 0; i < x.length; i++) {
				final String[] words = input.readLine().trim().split(" ");

				if (words[0].equalsIgnoreCase("pos")) {
					y[i] = 1;
				} else if (words[0].equalsIgnoreCase("neg")) {
					y[i] = 0;
				} else {
					System.err.println("Invalid class label: " + words[0]);
				}

				x[i] = words;
			}
		} catch (final IOException ex) {
			System.err.println(ex);
		}

		moviex = new double[x.length][];
		moviey = new int[y.length];
		final Bag<String> bag = new Bag<>(feature);
		for (int i = 0; i < x.length; i++) {
			moviex[i] = bag.feature(x[i]);
			moviey[i] = y[i];
		}
	}

	/**
	 * Test of predict method, of class NaiveBayes.
	 */
	@Test
	public void testPredict() {
		System.out.println("---predict---");
		final ArffParser arffParser = new ArffParser();
		arffParser.setResponseIndex(4);
		try {
			final AttributeDataset iris = arffParser.parse(new ClassPathResource("iris.arff").getFile());
			final double[][] x = iris.toArray(new double[iris.size()][]);
			final int[] y = iris.toArray(new int[iris.size()]);

			final int n = x.length;
			final LOOCV loocv = new LOOCV(n);
			int error = 0;
			for (int l = 0; l < n; l++) {
				final double[][] trainx = smile.math.Math.slice(x, loocv.train[l]);
				final int[] trainy = smile.math.Math.slice(y, loocv.train[l]);

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

						condprob[i][j] = new GaussianMixture(xi, 3);
					}
				}

				final NaiveBayes bayes = new NaiveBayes(priori, condprob);

				if (y[loocv.test[l]] != bayes.predict(x[loocv.test[l]]))
					error++;
			}

			System.out.format("Iris error rate = %.2f%%%n", 100.0 * error / x.length);
			assertEquals(8, error);
		} catch (final Exception ex) {
			System.err.println(ex);
		}
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnMultinomial() {
		System.out.println("---batch learn Multinomial---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.MULTINOMIAL, 2, feature.length);

			bayes.learn(trainx, trainy);

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);
			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Batch Multinomial error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnMultinomial2() {
		System.out.println("---online learn Multinomial---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.MULTINOMIAL, 2, feature.length);

			for (int j = 0; j < trainx.length; j++) {
				bayes.learn(trainx[j], trainy[j]);
			}

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);
			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Online Multinomial error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnPolyaUrn() {
		System.out.println("---batch learn PolyaUrn---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.POLYAURN, 2, feature.length);

			bayes.learn(trainx, trainy);

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);
			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Batch PolyaUrn error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnPolyaUrn2() {
		System.out.println("---online learn PolyaUrn---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.POLYAURN, 2, feature.length);

			for (int j = 0; j < trainx.length; j++) {
				bayes.learn(trainx[j], trainy[j]);
			}

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);
			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Online PolyaUrn error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 265);
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnBernoulli() {
		System.out.println("---batch learn Bernoulli---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.BERNOULLI, 2, feature.length);

			bayes.learn(trainx, trainy);

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);

			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Batch Bernoulli error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 270);
	}

	/**
	 * Test of learn method, of class SequenceNaiveBayes.
	 */
	@Test
	public void testLearnBernoulli2() {
		System.out.println("---online learn Bernoulli---");

		final double[][] x = moviex;
		final int[] y = moviey;
		final int n = x.length;
		final int k = 10;
		final CrossValidation cv = new CrossValidation(n, k);
		int error = 0;
		int total = 0;
		for (int i = 0; i < k; i++) {
			final double[][] trainx = smile.math.Math.slice(x, cv.train[i]);
			final int[] trainy = smile.math.Math.slice(y, cv.train[i]);
			final NaiveBayes bayes = new NaiveBayes(NaiveBayes.Model.BERNOULLI, 2, feature.length);

			for (int j = 0; j < trainx.length; j++) {
				bayes.learn(trainx[j], trainy[j]);
			}

			final double[][] testx = smile.math.Math.slice(x, cv.test[i]);
			final int[] testy = smile.math.Math.slice(y, cv.test[i]);

			for (int j = 0; j < testx.length; j++) {
				final int label = bayes.predict(testx[j]);
				if (label != -1) {
					total++;
					if (testy[j] != label) {
						error++;
					}
				}
			}
		}

		System.out.format("Online Bernoulli error = %d of %d, error rate = %.2f%% %n", error, total,
				100.0 * error / total);
		assertTrue(error < 270);
	}

}
