package edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian;

import java.util.List;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.log4j.Logger;
import org.jgap.Gene;
import org.jgap.IChromosome;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;

public class MultipleGaussianFit extends AbstractMultipleGaussianFit {
	private final static Logger log = Logger.getLogger(MultipleGaussianFit.class);
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;

	public MultipleGaussianFit(Double maxLimit) {
		super(maxLimit);
	}

	// http://en.wikipedia.org/wiki/Coefficient_of_determination
	@Override
	public double evaluate(IChromosome a_subject) {

		final Gene[] elements = a_subject.getGenes();
		gaussians.clear();
		for (int i = 0; i < elements.length; i = i + 3) {
			final MyGaussianFit gaussian = new MyGaussianFit((Double) elements[i].getAllele(),
					(Double) elements[i + 1].getAllele(), (Double) elements[i + 2].getAllele());
			gaussians.add(gaussian);
		}

		final double rSquare = evaluateFitting(experimentalData, gaussians, false);
		return rSquare;
	}

	public static double evaluateFitting(List<double[]> experimentalData, List<MyGaussianFit> gaussians,
			boolean printout) {
		// return calculateR2(experimentalData, gaussians, printout);
		return calculateRSS(experimentalData, gaussians, printout);
	}

	private static double calculateRSS(List<double[]> experimentalData, List<MyGaussianFit> fittedGaussians,
			boolean printout) {

		double RSS = 0.0;
		for (final double[] pair : experimentalData) {
			// if (pair[1] > 0) {
			final double expY = pair[1];
			double fittedY = 0.0;
			for (final MyGaussianFit gaussian : fittedGaussians) {
				fittedY += gaussian.getGaussianY(pair[0]);
			}
			RSS += Math.abs(expY - fittedY);
		}

		final double ret = 1 / RSS;
		return ret;
	}

	/**
	 * Calculates the R squared as followed in
	 * https://newonlinecourses.science.psu.edu/stat501/node/255/
	 * 
	 * @param experimentalData
	 * @param fittedGaussians
	 * @return
	 */

	private static double calculateR2(List<double[]> experimentalData, List<MyGaussianFit> fittedGaussians,
			boolean printout) {

		final SimpleRegression regression = new SimpleRegression();
		for (final double[] pair : experimentalData) {
			// if (pair[1] > 0) {
			final double expY = pair[1];
			double fittedY = 0.0;
			for (final MyGaussianFit gaussian : fittedGaussians) {
				fittedY += gaussian.getGaussianY(pair[0]);
			}
			regression.addData(expY, fittedY);
			if (printout) {
				System.out.println(expY + "\t" + fittedY);
			}
			// }
		}

		final double rSquare = regression.getRSquare() * 100;
		return rSquare;
	}

}
