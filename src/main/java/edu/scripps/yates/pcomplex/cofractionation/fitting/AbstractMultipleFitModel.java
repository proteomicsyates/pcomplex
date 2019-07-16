package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.util.ArrayList;
import java.util.List;

import org.jgap.Configuration;
import org.jgap.InvalidConfigurationException;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian.AbstractMultipleGaussianFit;
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
	private TDoubleList profile;
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

	public void setExperimentalData(TDoubleList profile) {
		this.profile = profile;

	}

	protected List<double[]> transformExperimentalDataForModelling(TDoubleList profile) {
		final List<double[]> ret = new ArrayList<double[]>();
		int x = 1;
		for (final double y : profile.toArray()) {
			final double[] point = new double[2];
			point[0] = x;
			point[1] = y;
			ret.add(point);
			x++;
		}
		return ret;
	}

	protected List<int[]> transformExperimentalDataForModelling(TIntList profile) {
		final List<int[]> ret = new ArrayList<int[]>();
		int x = 1;
		for (final int y : profile.toArray()) {
			final int[] point = new int[2];
			point[0] = x;
			point[1] = y;
			ret.add(point);
			x++;
		}
		return ret;
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
			experimentalData = transformExperimentalDataForModelling(profile);
		}
		return experimentalData;
	}

	public TDoubleList getProfile() {
		return profile;
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
