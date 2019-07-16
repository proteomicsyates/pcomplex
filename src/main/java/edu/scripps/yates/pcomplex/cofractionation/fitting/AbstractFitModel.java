package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.util.List;
import java.util.Set;

import org.jgap.Configuration;
import org.jgap.InvalidConfigurationException;
import org.jgap.impl.DefaultConfiguration;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian.AbstractGaussianFit;

public abstract class AbstractFitModel {
	protected static final double CHANGE_FITNESS_THRESHOLD = 0.0001;
	protected static final int MAX_ROUNDS_WITH_NO_CHANGE = 15;
	public static final String FITTING_ERROR = "fitting error";
	public static final double MINIMUM_VALID_RSQUARED = 0.5;
	List<double[]> experimentalData;
	protected final Configuration conf;
	protected AbstractGaussianFit optimizationFunction;
	protected MyGaussianFit fittedGaussian;
	protected static int progress = 0;
	protected Double fixedHeight1;
	protected Double fixedAvgBkg;
	protected Double fixedStdBkg;

	public AbstractFitModel(Double maxLimit) throws InvalidConfigurationException {
		Configuration.reset();
		conf = new DefaultConfiguration();

		optimizationFunction = createOptimizationFunction(maxLimit);
		conf.setFitnessFunction(optimizationFunction);

		conf.setPopulationSize(500);
	}

	public abstract AbstractGaussianFit createOptimizationFunction(Double maxLimit);

	public void setExperimentalData(List<double[]> optimization_set) {
		experimentalData = optimization_set;
	}

	public abstract Set<MyGaussianFit> runFit() throws InvalidConfigurationException;

	public abstract int getMaxIterations();

	/**
	 * @param fixedA the fixedA to set
	 */
	public void setFixedHeight(Double fixedA) {
		this.fixedHeight1 = fixedA;
	}

	/**
	 * @param fixedAvgBkg the fixedAvgBkg to set
	 */
	public void setFixedAverage(Double fixedAvgBkg) {
		this.fixedAvgBkg = fixedAvgBkg;
	}

	/**
	 * @param fixedStdBkg the fixedStdBkg to set
	 */
	public void setFixedStdev(Double fixedStdBkg) {
		this.fixedStdBkg = fixedStdBkg;
	}

	protected double getMaxY() {

		double max = -1;

		for (int i = 0; i < experimentalData.size(); i++) {
			if (experimentalData.get(i)[1] > max) {
				max = experimentalData.get(i)[1];
			}
		}

		return max;
	}

	public MyGaussianFit getFittedGaussian() {
		return fittedGaussian;
	}

}
