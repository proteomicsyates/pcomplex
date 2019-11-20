package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.util.List;

import org.jgap.Configuration;
import org.jgap.InvalidConfigurationException;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian.AbstractMultipleGaussianFit;
import edu.scripps.yates.pcomplex.util.PComplexUtil;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;

public abstract class AbstractMultipleFitModel {
	protected static final double CHANGE_FITNESS_THRESHOLD = 0.0001;
	protected static final int MAX_ROUNDS_WITH_NO_CHANGE = 15;
	public static final String FITTING_ERROR = "fitting error";
	public static final double MINIMUM_VALID_RSQUARED = 0.5;
	protected Configuration conf;
	protected AbstractMultipleGaussianFit optimizationFunction;
	protected List<MyGaussianFit> fittedGaussians;
	protected static int progress = 0;
	protected Double[] fixedHeight;
	protected Double[] fixedAvg;
	protected Double[] fixedStd;

	protected final int numGaussians;
	private TDoubleList processedProfile;
	private List<double[]> experimentalData;
	protected String acc;
	private TDoubleList rawProfile;
	private TIntList spcProfile;

	public AbstractMultipleFitModel(String acc, Double maxLimit, int numGaussians) {
		this.acc = acc;
		this.numGaussians = numGaussians;
		fixedHeight = new Double[numGaussians];
		fixedAvg = new Double[numGaussians];
		fixedStd = new Double[numGaussians];

		optimizationFunction = createOptimizationFunction(maxLimit);

	}

	public abstract AbstractMultipleGaussianFit createOptimizationFunction(Double maxLimit);

	public void setProcessedProfile(TDoubleList profile) {
		this.processedProfile = profile;

	}

	protected abstract ModelPlot runFit() throws InvalidConfigurationException;

	public abstract int getMaxIterations();

	/**
	 * @param fixedA
	 *            the fixedA to set
	 */
	public void setFixedHeight(Double fixedA, int gaussianIndex) {
		fixedHeight[gaussianIndex] = fixedA;
	}

	/**
	 * @param fixedAvgBkg
	 *            the fixedAvgBkg to set
	 */
	public void setFixedAverage(Double fixedAvgBkg, int gaussianIndex) {
		fixedAvg[gaussianIndex] = fixedAvgBkg;
	}

	/**
	 * @param fixedStdBkg
	 *            the fixedStdBkg to set
	 */
	public void setFixedStdev(Double fixedStdBkg, int gaussianIndex) {
		fixedStd[gaussianIndex] = fixedStdBkg;
	}

	protected double getMaxY() {

		double max = -1;
		final List<double[]> experimentalData = getExperimentalData();
		for (int i = 0; i < experimentalData.size(); i++) {
			if (experimentalData.get(i)[1] > max) {
				max = experimentalData.get(i)[1];
			}
		}

		return max;
	}

	public List<MyGaussianFit> getFittedGaussians() {
		return fittedGaussians;
	}

	protected List<double[]> getExperimentalData() {
		if (experimentalData == null) {
			experimentalData = PComplexUtil.transformExperimentalDataForModelling(processedProfile);
		}
		return experimentalData;
	}

	public TDoubleList getProcessedProfile() {
		return processedProfile;
	}

	public abstract double getFitnessValue();

	public TDoubleList getRawProfile() {
		return rawProfile;
	}

	public TIntList getSPCProfile() {
		return spcProfile;
	}

	public void setRawProfile(TDoubleList rawProfile) {
		this.rawProfile = rawProfile;
	}

	public void setSPCProfile(TIntList spcProfile) {
		this.spcProfile = spcProfile;
	}

	public abstract double getRSquare();
}
