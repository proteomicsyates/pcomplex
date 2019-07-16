package edu.scripps.yates.pcomplex;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.log4j.Logger;
import org.jgap.InvalidConfigurationException;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.FittingUtil;
import edu.scripps.yates.pcomplex.cofractionation.fitting.ModelPlot;
import edu.scripps.yates.pcomplex.cofractionation.fitting.MultipleGaussianFitModel;
import edu.scripps.yates.pcomplex.cofractionation.training.ClassificationResult;
import edu.scripps.yates.pcomplex.cofractionation.training.NaiveBayesSmile;
import edu.scripps.yates.pcomplex.cofractionation.training.NaiveBayesWeka;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.gui.ChartJFrame;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import smile.netlib.NLMatrix;

public class CoFractionationAnalyzer {
	private final static Logger log = Logger.getLogger(CoFractionationAnalyzer.class);

	public static void main(String[] args) {
		try {
			final File projectSummaryFile = new File(args[0]);
			final File folderForMatrixes = new File(args[1]);
			if (projectSummaryFile.isFile()) {

				final CoFractionationAnalyzer ca = new CoFractionationAnalyzer(projectSummaryFile, folderForMatrixes);
				ca.run();
			} else {
				final File[] listFiles = projectSummaryFile.listFiles();
				for (final File file : listFiles) {
					if (file.isFile()) {
						if (FilenameUtils.getExtension(file.getAbsolutePath()).equals("tsv")) {
							final CoFractionationAnalyzer ca = new CoFractionationAnalyzer(file, folderForMatrixes);
							ca.run();
						}
					}
				}
			}
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private final SeparationExperiment experiment;
	private final int minFractions = 5;
	private final int minSPCInOneFraction = 2;
	private final File logFile;
	private ProteinComplexDB referenceDB;
	private final File matrixFolder;
	private String projectName;
	private final boolean logFittings;
	private static final PearsonsCorrelation pearson = new PearsonsCorrelation();
	private static final EuclideanDistance euclideanDistance = new EuclideanDistance();
	private static final int PROFILE_SMOOTH_WIDTH = 3;
	private static final int MAX_NUM_GAUSSIANS = 4;
	// order to apply second filter on each region
	private static THashMap<String, List<MyGaussianFit>> fitsByProteins = new THashMap<String, List<MyGaussianFit>>();
	private static THashSet<String> profilesPrinted = new THashSet<String>();
	private final TLongArrayList fittingTimes = new TLongArrayList();

	private final ChartJFrame window = new ChartJFrame();

	private enum MLLIBRARY {
		SMILE, WEKA
	};

	private final MLLIBRARY machineLearningLibraryToUse = MLLIBRARY.WEKA;

	public CoFractionationAnalyzer(File projectSummaryFile, File folderForMatrixes) throws IOException {

		projectName = FilenameUtils.getBaseName(projectSummaryFile.getAbsolutePath());
		if (projectName.contains("_")) {
			projectName = projectName.substring(0, projectName.lastIndexOf("_"));
		}
		logFile = new File(projectSummaryFile.getParent() + File.separator + projectName + ".log");
		matrixFolder = folderForMatrixes;
		log.info("Reading data from project file");
		experiment = ProteinComplexAnalyzer.loadProjectSummaryFile(projectName, projectSummaryFile);
		log.info(experiment.getFractions().size() + " fractions in experiment " + projectName);
		logFittings = true;
	}

	public void run() throws Exception {
		log.info("Analyzing project " + projectName);
		final List<String> proteinList = new ArrayList<String>();

		final Map<String, List<Protein>> totalProteinsByAcc = experiment.getTotalProteinsByAcc(minFractions,
				minSPCInOneFraction);
		log.info(totalProteinsByAcc.size() + " proteins in total in all fractions");
		proteinList.addAll(totalProteinsByAcc.keySet());
		// it is important to sort the accessions, so that the order is kept
		Collections.sort(proteinList);
		// we are using NSAF values to build the elution profiles
		final DataType dataType = DataType.NSAF;

		final Map<DistanceMeasure, NLMatrix> distanceMap = new THashMap<DistanceMeasure, NLMatrix>();
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			distanceMap.put(distance, new NLMatrix(proteinList.size(), proteinList.size(), Double.NaN));
		}
		// check whether the matrixes files exist or not
		boolean calculateDistanceMeasurements = false;
		final List<Boolean> calculateDistance = new ArrayList<Boolean>(DistanceMeasure.values().length);
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			final File outputFileNameForMatrix = getOutputFileNameForMatrix(distance);
			if (!outputFileNameForMatrix.exists() || outputFileNameForMatrix.length() == 0l) {
				calculateDistanceMeasurements = true;
				calculateDistance.add(true);
			} else {
				calculateDistance.add(false);
			}
		}
		if (calculateDistanceMeasurements) {
			calculateDistanceMeasurements(proteinList, distanceMap, dataType, calculateDistance);
			// print files
			for (final DistanceMeasure distance : DistanceMeasure.values()) {
				printMatrix(proteinList, distanceMap.get(distance), distance);
			}
		} else {
			loadDistancesMatrixes(distanceMap);
		}
		preProcessMatrixes(distanceMap);
		final ClassificationResult result = naiveBayesTraining(distanceMap, proteinList);
		processResult(result);
	}

	private void processResult(ClassificationResult result) {
		// sort by the probability of being intra-complex and calculate the
		// TP/(TP+FP)
result.
	}

	private ClassificationResult naiveBayesTraining(Map<DistanceMeasure, NLMatrix> distanceMap,
			List<String> proteinList) throws Exception {
		if (machineLearningLibraryToUse == MLLIBRARY.SMILE) {
			final NaiveBayesSmile naiveBayesSmile = new NaiveBayesSmile() {

				@Override
				public ClassLabel getClassLabel(String protein, String protein2) throws IOException {
					return CoFractionationAnalyzer.this.getClassLabel(protein, protein2);
				}
			};
			return naiveBayesSmile.naiveBayesTraining(distanceMap, proteinList);
		} else {

			final String arfFileName = matrixFolder.getAbsolutePath() + File.separator + projectName + ".arff";
			final NaiveBayesWeka naiveBayesWeka = new NaiveBayesWeka(projectName, arfFileName) {

				@Override
				public ClassLabel getClassLabel(String protein, String protein2) throws IOException {
					return CoFractionationAnalyzer.this.getClassLabel(protein, protein2);
				}
			};
			return naiveBayesWeka.naiveBayesTraining(distanceMap, proteinList);
		}

	}

	private ClassLabel getClassLabel(String protein, String protein2) throws IOException {
		getReferenceDB();
		final Set<ProteinComplex> proteinComplexes1 = referenceDB.getProteinComplexesByProtein(protein);
		final Set<ProteinComplex> proteinComplexes2 = referenceDB.getProteinComplexesByProtein(protein2);
		if (proteinComplexes1.isEmpty() || proteinComplexes2.isEmpty()) {
			return ClassLabel.NOVEL_INTERACTORS;
		}
		for (final ProteinComplex proteinComplex : proteinComplexes1) {
			if (proteinComplexes2.contains(proteinComplex)) {
				return ClassLabel.INTRA_COMPLEX;
			}
		}
		return ClassLabel.INTER_COMPLEX;
	}

	private void loadDistancesMatrixes(Map<DistanceMeasure, NLMatrix> distanceMap) throws IOException {
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			final File matrixFile = getOutputFileNameForMatrix(distance);
			final NLMatrix matrix = loadDistanceMatrix(matrixFile);
			distanceMap.put(distance, matrix);
		}
		log.info(distanceMap.size() + " matrixes loaded");
	}

	private NLMatrix loadDistanceMatrix(File matrixFile) throws IOException {
		final long t1 = System.currentTimeMillis();
		log.info("Loading file with matrix of distances " + matrixFile.getAbsolutePath());
		final FileInputStream fstream = new FileInputStream(matrixFile);
		final BufferedReader br = new BufferedReader(new InputStreamReader(fstream));

		String strLine;
		NLMatrix matrix = null;
		// Read File Line By Line
		// each line is a new row in the matrix
		int row = 0;
		while ((strLine = br.readLine()) != null) {
			try {
				// skip line 0
				if (row == 0) {
					continue;
				}
				final String[] split = strLine.trim().split("\\s+");
				if (matrix == null) {
					final int matrixSize = split.length - 1;
					matrix = new NLMatrix(matrixSize, matrixSize, Double.NaN);
				}
				for (int col = 1; col < split.length; col++) {
					final double num = Double.valueOf(split[col]);
					matrix.set(row - 1, col - 1, num);
				}
			} finally {
				row++;
			}
		}

		// Close the input stream
		fstream.close();
		final long t2 = System.currentTimeMillis();
		log.info("Matrix of size " + matrix.nrows() + "x" + matrix.ncols() + " has been loaded in "
				+ DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));
		return matrix;
	}

	private void preProcessMatrixes(Map<DistanceMeasure, NLMatrix> distanceMap) {
		log.info("Preprocessing matrixes");
		// standardize values in matrixes
		// making mean=0 and standard deviation 1
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			standardizeMatrix(distanceMap.get(distance));
		}
		log.info(DistanceMeasure.values().length + " matrixes standardized");
	}

	private void calculateDistanceMeasurements(List<String> proteinList, Map<DistanceMeasure, NLMatrix> distanceMap,
			DataType dataType, List<Boolean> calculateDistance) throws IOException, InvalidConfigurationException {
		final FileWriter logFileWriter = new FileWriter(logFile, false);
		final long total = proteinList.size() * (proteinList.size() - 1l) / 2;
		final ProgressCounter counter = new ProgressCounter(total, ProgressPrintingType.PERCENTAGE_STEPS, 4);
		counter.setShowRemainingTime(true);
		for (int i = 0; i < proteinList.size(); i++) {
			final String protein1 = proteinList.get(i);
			for (int j = i + 1; j < proteinList.size(); j++) {
				final String protein2 = proteinList.get(j);
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info(printIfNecessary);
					final double averageTime = Maths.mean(fittingTimes);
					log.info(DatesUtil.getDescriptiveTimeFromMillisecs(averageTime) + " in average per protein");
				}
				if (calculateDistance.get(DistanceMeasure.CORRELATION.ordinal())) {
					final double correlationCoefficient = getPearsonCorrelationCoefficient(experiment, protein1,
							protein2, dataType);
					distanceMap.get(DistanceMeasure.CORRELATION).set(i, j, correlationCoefficient);
					log.debug("correlationCoefficient" + "\t" + correlationCoefficient);
				}
				if (calculateDistance.get(DistanceMeasure.CORRELATION_PVALUE.ordinal())) {
					final double pvalue = getCorrelationPValue(experiment, protein1, protein2, dataType);
					distanceMap.get(DistanceMeasure.CORRELATION_PVALUE).set(i, j, pvalue);
					log.debug("pvalue" + "\t" + pvalue);
				}
				if (calculateDistance.get(DistanceMeasure.EUCLIDEAN_DISTANCE.ordinal())) {
					final double euDistance = getEuclideanDistance(experiment, protein1, protein2, dataType);
					distanceMap.get(DistanceMeasure.EUCLIDEAN_DISTANCE).set(i, j, euDistance);
					log.debug("euDistance" + "\t" + euDistance);
				}
				if (calculateDistance.get(DistanceMeasure.PEAK_LOCATION.ordinal())) {
					final double peakLocation = getPeakLocation(experiment, protein1, protein2, dataType);
					log.debug("peakLocation" + "\t" + peakLocation);
					distanceMap.get(DistanceMeasure.PEAK_LOCATION).set(i, j, peakLocation);
				}
				if (calculateDistance.get(DistanceMeasure.COAPEX_SCORE.ordinal())) {
					final double coApexScore = getCoApexScore(experiment, protein1, protein2, dataType, logFileWriter,
							logFittings);
					if (i == 0 && j == 213) {
						log.info("asdf");
					}
					distanceMap.get(DistanceMeasure.COAPEX_SCORE).set(i, j, coApexScore);
					log.debug("coApexScore" + "\t" + coApexScore);
				}

			}
		}
	}

	private File getOutputFileNameForMatrix(DistanceMeasure distanceMeasure) {
		return new File(matrixFolder.getAbsolutePath() + File.separator + "Matrix_" + distanceMeasure.name() + "_"
				+ projectName + ".txt");
	}

	private void printMatrix(List<String> proteinList, NLMatrix matrix, DistanceMeasure distanceMeasure)
			throws IOException {
		log.info("Printing matrix for distance: " + distanceMeasure.name());
		final File fileOutput = getOutputFileNameForMatrix(distanceMeasure);
		final FileWriter fw = new FileWriter(fileOutput);

		for (final String acc : proteinList) {
			fw.write("\t" + acc);
		}
		fw.write("\n");
		for (int row = 0; row < matrix.nrows(); row++) {
			fw.write(proteinList.get(row) + "\t");
			for (int col = 0; col < matrix.ncols(); col++) {
				fw.write(matrix.get(row, col) + "\t");
			}
			fw.write("\n");
		}

		fw.close();
	}

	private ProteinComplexDB getReferenceDB() throws IOException {
		if (referenceDB == null) {
			ProteinComplexAnalyzer.useCoreCorumDB = true;
			ProteinComplexAnalyzer.useComplexPortalDB = false;
			ProteinComplexAnalyzer.useHUMAP = false;
			referenceDB = ProteinComplexAnalyzer.getDBs().get(0);
		}
		return referenceDB;
	}

	/**
	 * The difference, in fractions, between the locations of the maximum values
	 * of the fractionation profiles of the two proteins
	 * 
	 * @param exp
	 * @param acc1
	 * @param acc2
	 * @param dataType
	 * @return
	 */
	public int getPeakLocation(SeparationExperiment exp, String acc1, String acc2, DataType dataType) {
		final TDoubleArrayList profile1 = getDataProfile(exp, acc1, dataType);
		final TDoubleArrayList profile2 = getDataProfile(exp, acc2, dataType);

		// get the number of the fractions in which there is the maximum peak of
		// the
		// profile.
		// it is a list just in case we have the maximum in multiple fractions
		final TIntList maximums1 = getFractionsWithMaximums(profile1);
		final TIntList maximums2 = getFractionsWithMaximums(profile2);

		// now, get the minimum distance between both
		final int distance = getMinimumDistance(maximums1, maximums2);
		return distance;
	}

	private int getMinimumDistance(TIntList array1, TIntList array2) {
		int ret = Integer.MAX_VALUE;
		for (final int fraction1 : array1.toArray()) {
			for (final int fraction2 : array2.toArray()) {
				final int abs = Math.abs(fraction1 - fraction2);
				if (abs < ret) {
					ret = abs;
				}
			}
		}
		return ret;
	}

	private void printFitting(List<MyGaussianFit> fits, TDoubleArrayList profile, String acc, FileWriter fw)
			throws IOException {
		if (profilesPrinted.contains(acc)) {
			return;
		}
		profilesPrinted.add(acc);
		fw.write("\nFitting for " + acc + "\n");
		for (int fractionNumber = 0; fractionNumber < profile.size(); fractionNumber++) {
			fw.write(fractionNumber + "\t");
		}
		fw.write("\n");
		for (final double data : profile.toArray()) {
			fw.write(data + "\t");
		}
		fw.write("\n");
		for (final MyGaussianFit fitting : fits) {
			for (int fractionNumber = 0; fractionNumber < profile.size(); fractionNumber++) {
				fw.write(fitting.getGaussian().value(fractionNumber) + "\t");
			}
			fw.write("\n");
		}
		fw.flush();
	}

	private List<MyGaussianFit> fitToGaussiansWithGeneticAlgorithm(TDoubleArrayList profile, String acc,
			TIntList spcList) throws InvalidConfigurationException {
		if (fitsByProteins.containsKey(acc)) {
			return fitsByProteins.get(acc);
		}
		final long t0 = System.currentTimeMillis();
		// fill missing values with the average of neighbour points
		final TDoubleList profileWithImputedMissingValues = imputeMissingValuesOfProfile(profile);

		// smooth profile by a sliding average of width 5
		final TDoubleList profileSmoothed = smoothProfile(profileWithImputedMissingValues, PROFILE_SMOOTH_WIDTH);

		log.info("Fitting profile for protein " + acc);
		if (acc.equals("A5A3E0")) {
			log.info(acc);
		}
		// here I will store the list of fits, with the gaussians and the error
		// of the
		// fit
		double maxFitnessValue = -Double.MAX_VALUE;
		final List<MyGaussianFit> ret = new ArrayList<MyGaussianFit>();
		for (int numGaussians = 1; numGaussians <= MAX_NUM_GAUSSIANS; numGaussians++) {
			final ModelPlot modelPlot = launchFitting(acc, numGaussians, spcList, profile, profileSmoothed);
			final double fitnessValue = modelPlot.getFitnessValue();

			if (fitnessValue > maxFitnessValue) {
				FittingUtil.printFitness(System.out, modelPlot.getFitModel().getFittedGaussians(),
						modelPlot.getFitnessValue());
				Double previousFitnessValue = maxFitnessValue;
				if (previousFitnessValue == -Double.MAX_VALUE) {
					previousFitnessValue = null;
				}
				window.showGraph(modelPlot, previousFitnessValue);
				maxFitnessValue = fitnessValue;
				ret.clear();
				ret.addAll(modelPlot.getFitModel().getFittedGaussians());
			}
		}
		fitsByProteins.put(acc, ret);

		final long t1 = System.currentTimeMillis();
		fittingTimes.add(t1 - t0);
		return ret;
	}

	private ModelPlot launchFitting(String acc, int numGaussians, TIntList spcProfile, TDoubleList rawProfile,
			TDoubleList profile) throws InvalidConfigurationException {
		final MultipleGaussianFitModel model = new MultipleGaussianFitModel(acc, numGaussians);

		model.setExperimentalData(profile);
		model.setRawProfile(rawProfile);
		model.setSPCProfile(spcProfile);
		return model.runFit();
	}

	private TDoubleList smoothProfile(TDoubleList profile, int slideWindow) {
		final int offsetMin = (slideWindow - 1) / 2;
		final TDoubleList ret = new TDoubleArrayList(profile.size());
		for (int index = 0; index < profile.size(); index++) {
			final TDoubleArrayList toAverage = new TDoubleArrayList();
			for (int offset = -offsetMin; offset <= offsetMin; offset++) {
				if (index + offset >= 0 && index + offset < profile.size()) {
					toAverage.add(profile.get(index + offset));
				}
			}
			ret.add(Maths.mean(toAverage));
		}
		return ret;
	}

	/**
	 * It imputes the missing values of the profile by averaging them with their
	 * neighbors points
	 * 
	 * @param profile
	 * @return
	 */
	private TDoubleList imputeMissingValuesOfProfile(TDoubleArrayList profile) {
		final TDoubleList ret = new TDoubleArrayList(profile.size());
		for (int index = 0; index < profile.size(); index++) {
			// if (index == 38) {
			// log.info("asdf");
			// }
			final double point = profile.get(index);
			if (Double.compare(0.0, point) == 0) {
				double previousPoint = 0.0;
				double nextPoint = 0.0;
				if (index - 1 >= 0) {
					previousPoint = profile.get(index - 1);
				}
				if (index + 1 < profile.size()) {
					nextPoint = profile.get(index + 1);
				}
				final double newPoint = (previousPoint + nextPoint) / 2.0;
				ret.add(newPoint);
			} else {
				ret.add(point);
			}
		}
		return ret;
	}

	private void printPoints(List<WeightedObservedPoint> obs) {
		for (final WeightedObservedPoint weightedObservedPoint : obs) {
			System.out.print(weightedObservedPoint.getY() + ", ");
		}
		System.out.println();
	}

	/**
	 * Returns the sum of the squares of the distances between the observed and
	 * the calculated points
	 * 
	 * @param gaussian
	 * @param obs
	 * @return
	 */
	private double getError(Gaussian gaussian, List<WeightedObservedPoint> obs) {
		double error = 0.0;
		for (final WeightedObservedPoint point : obs) {
			final double x = point.getX();
			final double calculatedY = gaussian.value(x);
			final double observedY = point.getY();
			final double pointError = calculatedY - observedY;
			error += Math.pow(pointError, 2);
		}
		return error;
	}

	/**
	 * Returns the sum of the squares of the distances between the entire
	 * profile and the calculated points, but only taking into account data
	 * entering on the range of +- 3 sigmas (99.7% of the data)
	 * 
	 * @param gaussian
	 * @param obs
	 * @return
	 */
	private double getError(Gaussian gaussian, double mean, double sigma, TDoubleList obs) {
		double error = 0.0;
		final double rangeMin = mean - sigma * 3.0;
		final double rangeMax = mean + sigma * 3.0;
		for (int numFraction = 0; numFraction < obs.size(); numFraction++) {
			final double calculatedY = gaussian.value(numFraction);
			if (numFraction < rangeMin || numFraction > rangeMax) {
				continue;
			}
			final double observedY = obs.get(numFraction);

			final double pointError = calculatedY - observedY;
			error += Math.pow(pointError, 2);
		}
		return error;
	}

	/**
	 * Take the two gaussian fits that are closer in the profiles (closest mean)
	 * and calculate the euclidean distance between (mu,sigma) vs (mu,sigma) of
	 * the 2 gaussians.
	 * 
	 * @param fits1
	 *            an array result of fitting a Gaussian, with first element as
	 *            height, second as mean and third as sigma
	 * @param fits2
	 *            same as fits2 for the second Gaussian
	 * @return
	 */
	private double getDistanceBetweenClosestGaussians(List<MyGaussianFit> fits1, List<MyGaussianFit> fits2) {
		Pair<TDoubleList, TDoubleList> closestGaussians = null;
		double minDistance = Double.MAX_VALUE;
		for (final MyGaussianFit gaussian1 : fits1) {
			for (final MyGaussianFit gaussian2 : fits2) {
				final double mean1 = FittingUtil.getMean(gaussian1.getFittedParameters());// mean
				final double mean2 = FittingUtil.getMean(gaussian2.getFittedParameters());// mean
				final double distance = Math.abs(mean1 - mean2);
				if (distance < minDistance) {
					minDistance = distance;
					closestGaussians = new Pair<TDoubleList, TDoubleList>(gaussian1.getFittedParameters(),
							gaussian2.getFittedParameters());
				}
			}
		}
		if (closestGaussians == null) {
			log.info(closestGaussians);
		}
		final double distance = getEuclideanDistances(closestGaussians);
		return distance;
	}

	private double getEuclideanDistances(Pair<TDoubleList, TDoubleList> closestGaussians) {
		final TDoubleList array1 = new TDoubleArrayList();
		array1.add(FittingUtil.getMean(closestGaussians.getFirstelement()));// mean
		array1.add(FittingUtil.getSigma(closestGaussians.getFirstelement()));// sigma

		final TDoubleList array2 = new TDoubleArrayList();
		array2.add(FittingUtil.getMean(closestGaussians.getSecondElement()));// mean
		array2.add(FittingUtil.getSigma(closestGaussians.getSecondElement()));// sigma

		final double distance = euclideanDistance.compute(array1.toArray(), array2.toArray());
		return distance;
	}

	/**
	 * get the number of the fractions in which there is the maximum peak of the
	 * profile. <br>
	 * It is a list just in case we have the maximum in multiple fractions
	 * 
	 * @param profile1
	 * @return
	 */
	private TIntList getFractionsWithMaximums(TDoubleArrayList profile1) {
		final TIntList ret = new TIntArrayList();
		final double max = profile1.max();

		for (int fraction = 0; fraction < profile1.size(); fraction++) {
			if (Double.compare(profile1.get(fraction), max) == 0) {
				ret.add(fraction);
			}
		}
		return ret;
	}

	public double getPearsonCorrelationCoefficient(SeparationExperiment exp, String acc1, String acc2,
			DataType dataType) {
		final TDoubleArrayList profile1 = getDataProfile(exp, acc1, dataType);
		final TDoubleArrayList profile2 = getDataProfile(exp, acc2, dataType);

		final double correlation = pearson.correlation(profile1.toArray(), profile2.toArray());
		log.debug("Pearson correlation between " + acc1 + " - " + acc2 + " =\t" + correlation);

		return 1 - correlation;
	}

	/**
	 * Euclidean distance between the closest (mu,sigma) pairs, where mu and
	 * sigma are Gaussian parameters fitted to both fractionation profiles
	 * 
	 * @param exp
	 * @param acc1
	 * @param acc2
	 * @param dataType
	 * @return
	 * @throws IOException
	 * @throws InvalidConfigurationException
	 */
	public double getCoApexScore(SeparationExperiment exp, String acc1, String acc2, DataType dataType, FileWriter fw,
			boolean logFittings) throws IOException, InvalidConfigurationException {
		final TDoubleArrayList profile1 = getDataProfile(exp, acc1, dataType);
		final TIntList spcprofile1 = getDataProfileIntegers(exp, acc1, DataType.SPC);
		final TDoubleArrayList profile2 = getDataProfile(exp, acc2, dataType);
		final TIntList spcprofile2 = getDataProfileIntegers(exp, acc2, DataType.SPC);
		final List<MyGaussianFit> fits1 = fitToGaussiansWithGeneticAlgorithm(profile1, acc1, spcprofile1);
		if (logFittings) {
			printFitting(fits1, profile1, acc1, fw);
		}
		if (fits1.isEmpty()) {
			throw new IllegalArgumentException("This shoudn't happen");
		}

		final List<MyGaussianFit> fits2 = fitToGaussiansWithGeneticAlgorithm(profile2, acc2, spcprofile2);
		if (logFittings) {
			printFitting(fits2, profile2, acc2, fw);
		}
		if (fits2.isEmpty()) {
			throw new IllegalArgumentException("This shoudn't happen");
		}
		final double minimumDistance = getDistanceBetweenClosestGaussians(fits1, fits2);
		return minimumDistance;
	}

	private TDoubleArrayList getDataProfile(SeparationExperiment exp, String acc1, DataType dataType) {
		final TDoubleArrayList ret = new TDoubleArrayList();
		for (final Fraction fraction : exp.getSortedFractions()) {
			final Protein protein = fraction.getProteinByAcc(acc1);
			if (protein != null) {
				ret.add(protein.getData(dataType));
			} else {
				ret.add(0.0);
			}
		}
		return ret;
	}

	private TIntList getDataProfileIntegers(SeparationExperiment exp, String acc1, DataType dataType) {
		final TIntList ret = new TIntArrayList();
		for (final Fraction fraction : exp.getSortedFractions()) {
			final Protein protein = fraction.getProteinByAcc(acc1);
			if (protein != null) {
				ret.add(Double.valueOf(protein.getData(dataType)).intValue());
			} else {
				ret.add(0);
			}
		}
		return ret;
	}

	public double getCorrelationPValue(SeparationExperiment exp, String acc1, String acc2, DataType dataType) {
		final double[][] rectangularData = getRectangularData(exp, acc1, acc2, dataType);
		final PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation(rectangularData);
		final RealMatrix matrix = pearsonsCorrelation.getCorrelationPValues();
		final double entry = matrix.getEntry(0, 1);
		return entry;
	}

	public double getEuclideanDistance(SeparationExperiment exp, String acc1, String acc2, DataType dataType) {
		final TDoubleArrayList profile1 = getDataProfile(exp, acc1, dataType);
		final TDoubleArrayList profile2 = getDataProfile(exp, acc2, dataType);
		final double euDistance = euclideanDistance.compute(profile1.toArray(), profile2.toArray());
		return euDistance;
	}

	/**
	 * Get a matrix, which is two columns, one per each profile of each protein
	 * 
	 * @param exp
	 * @param acc1
	 * @param acc2
	 * @param dataType
	 * @return
	 */
	private double[][] getRectangularData(SeparationExperiment exp, String acc1, String acc2, DataType dataType) {
		final TDoubleArrayList profile1 = getDataProfile(exp, acc1, dataType);
		final TDoubleArrayList profile2 = getDataProfile(exp, acc2, dataType);
		final double[][] ret = new double[profile1.size()][2];
		for (int i = 0; i < profile1.size(); i++) {
			final double d1 = profile1.get(i);
			final double d2 = profile2.get(i);
			ret[i][0] = d1;
			ret[i][1] = d2;
		}
		return ret;
	}

	/**
	 * Take all non Nan values from the matrix and normalize them so that the
	 * mean is 0 and the standard deviation is 1
	 * 
	 * @param nlMatrix
	 */
	public void standardizeMatrix(NLMatrix matrix) {
		final long t1 = System.currentTimeMillis();
		log.info("Standardizing matrix of " + matrix.nrows() + "x" + matrix.ncols());
		final int cols = matrix.ncols();
		final int rows = matrix.nrows();
		final TDoubleArrayList numbers = new TDoubleArrayList();
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				final double num = matrix.get(row, col);
				if (!Double.isNaN(num)) {
					numbers.add(num);
				}
			}
		}
		final double mean = Maths.mean(numbers);
		final double sigma = Maths.stddev(numbers);
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				final double num = matrix.get(row, col);
				if (!Double.isNaN(num)) {
					final double normalizedNum = (num - mean) / sigma;
					matrix.set(row, col, normalizedNum);
				}
			}
		}
		final long t2 = System.currentTimeMillis();
		log.info("Matrix of " + matrix.nrows() + "x" + matrix.ncols() + " is now standard (mean=0,sigma=1). done in "
				+ DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));
	}

}
