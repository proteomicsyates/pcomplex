package edu.scripps.yates.pcomplex.cofractionation;

import org.apache.commons.math3.analysis.function.Gaussian;

import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

public class MyGaussianFit {
	private final Gaussian gaussian;
	private final TDoubleList fit;
	private double error;

	public MyGaussianFit(Gaussian gaussian, TDoubleList fit, double error) {
		this.gaussian = gaussian;
		this.fit = fit;
		this.error = error;
	}

	public MyGaussianFit(double h, double mean, double sigma) {
		fit = new TDoubleArrayList();
		fit.add(h);
		fit.add(mean);
		fit.add(sigma);
		gaussian = new Gaussian(h, mean, sigma);
	}

	public Gaussian getGaussian() {
		return gaussian;
	}

	public TDoubleList getFittedParameters() {
		return fit;
	}

	public double getError() {
		return error;
	}

	public double getGaussianY(double x) {
		return this.gaussian.value(x);
	}

}
