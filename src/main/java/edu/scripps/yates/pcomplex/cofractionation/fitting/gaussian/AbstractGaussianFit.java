package edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian;

import java.util.List;

import org.jgap.FitnessFunction;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;

public abstract class AbstractGaussianFit extends FitnessFunction {
	protected MyGaussianFit gaussian1;

	protected AbstractGaussianFit() {
	}

	/**
	 *
	 */
	private static final long serialVersionUID = -1893079012645806849L;
	protected List<double[]> experimentalData;

	public void setExperimentalData(List<double[]> experimentalData) {
		this.experimentalData = experimentalData;
	}

	protected abstract double evaluateFunction(double x);

}
