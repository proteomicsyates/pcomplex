package edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian;

import java.util.ArrayList;
import java.util.List;

import org.jgap.FitnessFunction;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;

public abstract class AbstractMultipleGaussianFit extends FitnessFunction {
	protected List<MyGaussianFit> gaussians = new ArrayList<MyGaussianFit>();
	protected final Double maxLimit;

	protected AbstractMultipleGaussianFit(Double maxLimit) {
		this.maxLimit = maxLimit;
	}

	/**
	 *
	 */
	private static final long serialVersionUID = -1893079012645806849L;
	protected List<double[]> experimentalData;

	public void setExperimentalData(List<double[]> experimentalData) {
		this.experimentalData = experimentalData;
	}

}
