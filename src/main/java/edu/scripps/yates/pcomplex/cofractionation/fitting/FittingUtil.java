package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jgap.Gene;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import gnu.trove.list.TDoubleList;

public class FittingUtil {

	public static final int MAX_GAUSSIAN_WIDTH = 10;// we consider this as the
													// maximum allowed width
													// (number of
													// fractions) in which a
													// protein in a complex can
													// elute;

	public static List<MyGaussianFit> getGaussiansFromGenes(Gene[] genes) {
		final List<MyGaussianFit> gaussians = new ArrayList<MyGaussianFit>();
		for (int i = 0; i < genes.length; i = i + 3) {
			final MyGaussianFit gaussian = new MyGaussianFit((Double) genes[i].getAllele(),
					(Double) genes[i + 1].getAllele(), (Double) genes[i + 2].getAllele());
			gaussians.add(gaussian);
		}
		return gaussians;
	}

	public static double getMean(TDoubleList fittedParameters) {
		return fittedParameters.get(1);
	}

	public static double getSigma(TDoubleList fittedParameters) {
		return fittedParameters.get(2);
	}

	public static double getHeight(TDoubleList fittedParameters) {
		return fittedParameters.get(0);
	}

	/**
	 * Returns mean + sigma*factor
	 * 
	 * @param fittedParameters
	 * @param factor
	 * @return
	 */
	public static double getXSigma(TDoubleList fittedParameters, double factor) {
		return getMean(fittedParameters) + getSigma(fittedParameters) * factor;
	}

	public static void printFitness(PrintStream out, List<MyGaussianFit> fittedGaussians, double rSquared2) {
		out.print("\nFitting with " + fittedGaussians.size() + " gaussians. R^2=" + rSquared2 + "\n");
		for (final MyGaussianFit myGaussianFit : fittedGaussians) {
			printGaussian(out, myGaussianFit);
		}
	}

	public static void printGaussian(PrintStream out, MyGaussianFit gaussian) {
		out.print("H=" + gaussian.getFittedParameters().get(0) + "\tMean=" + gaussian.getFittedParameters().get(1)
				+ "\tSigma=" + gaussian.getFittedParameters().get(2) + "\n");
	}

	public static double calculateR2(List<MyGaussianFit> fittedGaussians, List<double[]> experimentalData) {
		final SimpleRegression regression = new SimpleRegression();
		for (final double[] ds : experimentalData) {
			final double experimentalX = ds[0];
			final double experimentalY = ds[1];
			double gaussianY = 0.0;
			for (final MyGaussianFit gaussian : fittedGaussians) {
				gaussianY += gaussian.getGaussianY(experimentalX);
			}
			regression.addData(experimentalY, gaussianY);
		}
		final double rSquare = regression.getRSquare();
		return rSquare;
	}
}
