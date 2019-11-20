package edu.scripps.yates.pcomplex;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.log4j.Logger;
import org.jgap.InvalidConfigurationException;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.FittingCache;
import edu.scripps.yates.pcomplex.cofractionation.fitting.FittingUtil;
import edu.scripps.yates.pcomplex.cofractionation.fitting.ModelPlot;
import edu.scripps.yates.pcomplex.cofractionation.fitting.MultipleGaussianFitModel;
import edu.scripps.yates.pcomplex.cofractionation.training.ClassificationResult;
import edu.scripps.yates.pcomplex.cofractionation.training.CoFractionationDataset;
import edu.scripps.yates.pcomplex.cofractionation.training.MyClassifier;
import edu.scripps.yates.pcomplex.cofractionation.training.NaiveBayesWeka;
import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;
import edu.scripps.yates.pcomplex.cofractionation.training.TrueClassifier;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.distances.ApexScore;
import edu.scripps.yates.pcomplex.distances.EuclideanDistanceCalculator;
import edu.scripps.yates.pcomplex.distances.JaccardScore;
import edu.scripps.yates.pcomplex.distances.MutualInformation;
import edu.scripps.yates.pcomplex.distances.PearsonCorrelation;
import edu.scripps.yates.pcomplex.distances.PearsonCorrelationPlusNoise;
import edu.scripps.yates.pcomplex.distances.WeightedCrossCorrelation;
import edu.scripps.yates.pcomplex.gui.FittingProcessWindow;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.ClassificationEvaluation;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import edu.scripps.yates.pcomplex.util.DataType;
import edu.scripps.yates.pcomplex.util.ImageGenerator;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import smile.netlib.NLMatrix;

public class CoFractionationAnalyzer implements TrueClassifier {
	private final static Logger log = Logger.getLogger(CoFractionationAnalyzer.class);

	public static void main(String[] args) {
		try {
			final File projectSummaryFile = new File(args[0]);
			final File folderForMatrixes = new File(args[1]);
			final double precisionCutoff = 0.5;
			if (projectSummaryFile.isFile()) {

				final CoFractionationAnalyzer ca = new CoFractionationAnalyzer(projectSummaryFile, folderForMatrixes,
						precisionCutoff);
				ca.run();
			} else {
				final File[] listFiles = projectSummaryFile.listFiles();
				for (final File file : listFiles) {
					if (file.isFile()) {
						if (FilenameUtils.getExtension(file.getAbsolutePath()).equals("tsv")) {
							final CoFractionationAnalyzer ca = new CoFractionationAnalyzer(file, folderForMatrixes,
									precisionCutoff);
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
	private ProteinComplexDB referenceDBSimplified;
	private ProteinComplexDB referenceDB;
	private final File matrixFolder;
	private final String projectName;
	private final boolean logFittings;
	private static final int PROFILE_SMOOTH_WIDTH = 3;
	private static final boolean averageMissingValues = true;;

	private static final int MAX_NUM_GAUSSIANS = 4;
	public static final double SMALLEST_GAUSSIAN_PCT = 0.1;
	private static final int MAX_NUM_ITERATIONS_WITH_NO_IMPROVEMENT = 10;
	private static final double MIN_OVERLAPPING = 0.25; // overlapping to
														// consider 2 complexes
														// the same

	/*******
	 * REFERENCE SET PARAMETERS FOR LEARNING
	 */
	private static final double maxOverlapScoreInReferenceSet = 0.8;
	private static final int maxComplexSizeInReferenceSetForLearning = 50;
	private static final int minComplexSizeInReferenceSetForLearning = 3;
	private static final int kFoldCrossValidation = 2;
	private static final double minCorrelationForTraining = 0.5;
	/******************************/

	private static THashMap<String, List<MyGaussianFit>> fitsByProteins = new THashMap<String, List<MyGaussianFit>>();
	private static THashSet<String> profilesPrinted = new THashSet<String>();
	private final TLongArrayList fittingTimes = new TLongArrayList();

	private final FittingProcessWindow window = new FittingProcessWindow();

	private enum MLLIBRARY {
		SMILE, WEKA
	};

	private final MLLIBRARY machineLearningLibraryToUse = MLLIBRARY.WEKA;
	private final double procesionCutOff;
	private final boolean performMatchineLearning = true;
	private final boolean optimizeClusterOneMinDensity = false;

	public CoFractionationAnalyzer(File projectSummaryFile, File folderForMatrixes, double precisionCutoff)
			throws IOException {
		procesionCutOff = precisionCutoff;
		projectName = FilenameUtils.getBaseName(projectSummaryFile.getAbsolutePath());
		// if (projectName.contains("_")) {
		// projectName = projectName.substring(0, projectName.lastIndexOf("_"));
		// }
		logFile = new File(projectSummaryFile.getParent() + File.separator + projectName + ".log");
		matrixFolder = new File(folderForMatrixes.getAbsolutePath() + File.separator + getNameFromParams());
		log.info("Reading data from project file " + projectName);
		experiment = ProteinComplexAnalyzer.loadProjectSummaryFileNEW(projectName, projectSummaryFile);
		log.info(experiment.getFractions().size() + " fractions in experiment " + projectName);
		logFittings = true;
	}

	private String getNameFromParams() {
		final StringBuilder sb = new StringBuilder();
		// reference DB used
		if (ProteinComplexAnalyzer.useHUMAP) {
			StringUtils.addIfNotEmpty(sb, "_");
			sb.append("HUMAP");
		}
		if (ProteinComplexAnalyzer.useComplexPortalDB) {
			StringUtils.addIfNotEmpty(sb, "_");
			sb.append("COMPLEXPORTAL");
		}
		if (ProteinComplexAnalyzer.useCoreCorumDB) {
			StringUtils.addIfNotEmpty(sb, "_");
			sb.append("CORUM");
		}
		if (ProteinComplexAnalyzer.useAPID) {
			StringUtils.addIfNotEmpty(sb, "_");

			sb.append("APID");
		}
		//
		// parameters from reference database processing
		// maxOverlapScoreInReferenceSet
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("maxRefOv" + maxOverlapScoreInReferenceSet);
		// minComplexSizeInReferenceSetForLearning
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("minRefSiz" + minComplexSizeInReferenceSetForLearning);
		// maxComplexSizeInReferenceSetForLearning
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("maxRefSiz" + maxComplexSizeInReferenceSetForLearning);

		//
		// parameters for signal processing
		// PROFILE_SMOOTH_WIDTH
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("smooWd" + PROFILE_SMOOTH_WIDTH);

		// averageMissingValues
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("avgMisVal" + averageMissingValues);

		// minCorrelationForTraining
		StringUtils.addIfNotEmpty(sb, "_");
		sb.append("minCorTr" + minCorrelationForTraining);

		return sb.toString();
	}

	public void run() throws Exception {
		log.info("Analyzing project " + projectName);
		final List<String> proteinKeyList = new ArrayList<String>();

		final Map<String, List<Protein>> totalProteinsByAcc = experiment.getTotalProteinsByAcc(minFractions,
				minSPCInOneFraction);
		log.info(totalProteinsByAcc.size() + " proteins in total in all fractions");
		proteinKeyList.addAll(totalProteinsByAcc.keySet());
		// it is important to sort the accessions, so that the order is kept
		Collections.sort(proteinKeyList);
		// we are using NSAF values to build the elution profiles
		final DataType dataType = DataType.SPC;

		final Map<DistanceMeasure, NLMatrix> matrixesByDistanceMap = new THashMap<DistanceMeasure, NLMatrix>();
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			matrixesByDistanceMap.put(distance, new NLMatrix(proteinKeyList.size(), proteinKeyList.size(), Double.NaN));
		}
		// check whether the matrixes files exist or not
		boolean calculateDistanceMeasurements = false;
		final List<Boolean> calculateDistance = new ArrayList<Boolean>(DistanceMeasure.values().length);
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			// if (distance == DistanceMeasure.COAPEX_SCORE) {
			// calculateDistance.add(false);
			// continue;
			// }
			final File outputFileNameForMatrix = getOutputFileNameForMatrix(distance);
			if (!outputFileNameForMatrix.exists() || outputFileNameForMatrix.length() == 0l) {
				calculateDistanceMeasurements = true;
				calculateDistance.add(true);
			} else {
				calculateDistance.add(false);
			}
		}
		if (calculateDistanceMeasurements) {
			calculateDistanceMeasurements(proteinKeyList, matrixesByDistanceMap, dataType, calculateDistance);
			// print files
			for (final DistanceMeasure distance : DistanceMeasure.values()) {
				printMatrix(proteinKeyList, matrixesByDistanceMap.get(distance), distance);
			}
		} else {
			loadDistancesMatrixes(matrixesByDistanceMap);
		}
		if (!performMatchineLearning) {
			return;
		}

		// create dataset
		final CoFractionationDataset dataset = new CoFractionationDataset(matrixesByDistanceMap, proteinKeyList, this,
				minCorrelationForTraining, getReferenceDBSimplified(), matrixFolder, projectName);
		// pre process matrixes: make logs of probabilities, scale to have mean
		// 0 and sigma 1
		// preProcessMatrixes(distanceMap);
		// filter out proteins that are not in the reference database

		// train the machine learning
		final ClassificationResult result = naiveBayesTraining(dataset);
		createCurvesImages(result);
		final double pValueCutOff = estimatePValueThreshold(result, procesionCutOff);

		final List<ProteinComplex> finalClusters = clusterInteractions(result, pValueCutOff);
		if (finalClusters != null && !finalClusters.isEmpty()) {
			evaluateNovelyOfClusters(matrixFolder, finalClusters, getReferenceDB());
		} else {
			log.info("No clusters derived from experiment " + experiment.getProjectName());
		}
	}

	/**
	 * This function will say how many {@link ProteinComplex} are new and how
	 * many are already in the referenceDB
	 * 
	 * @param predictedComplexes
	 * @param referenceDB2
	 * @throws IOException
	 */
	private void evaluateNovelyOfClusters(File folder, List<ProteinComplex> predictedComplexes,
			ProteinComplexDB referenceDB2) throws IOException {
		final File fileGenes = new File(
				folder.getAbsolutePath() + File.separator + projectName + "_" + "complexes_GENES.txt");
		final File fileProteins = new File(
				folder.getAbsolutePath() + File.separator + projectName + "_" + "complexes_PROTEINS.txt");
		final FileWriter fwGenes = new FileWriter(fileGenes);
		final FileWriter fwProteins = new FileWriter(fileProteins);
		final Set<ProteinComplex> referenceComplexes = referenceDB2.getProteinComplexes();
		final Set<ProteinComplex> referenceComplexesSet = new THashSet<ProteinComplex>();
		referenceComplexesSet.addAll(referenceComplexes);
		fwGenes.write("Type" + "\t" + "Significance" + "\t" + "Density" + "\t" + "Max overlap with Ref" + "\t"
				+ "Reference Name" + "\t" + "Components" + "\n");
		for (final ProteinComplex predictedComplex : predictedComplexes) {
			final StringBuilder sbGenes = new StringBuilder();
			// calculate overlap
			double maxOverlap = 0.0;
			ProteinComplex bestReferenceComplex = null;
			for (final ProteinComplex referenceComplex : referenceComplexesSet) {
				final double overlap = ClusterEvaluation.getOverlap(predictedComplex, referenceComplex);
				if (overlap > maxOverlap) {
					maxOverlap = overlap;
					bestReferenceComplex = referenceComplex;
				}
			}
			// type
			if (maxOverlap > MIN_OVERLAPPING) {
				sbGenes.append("Known\t");
			} else {
				sbGenes.append("New\t");
			}
			// significance
			sbGenes.append(predictedComplex.getSignificance() + "\t");
			// density
			sbGenes.append(predictedComplex.getDensity() + "\t");
			// overlap
			sbGenes.append(maxOverlap + "\t");
			// best reference complex name
			if (bestReferenceComplex != null) {
				sbGenes.append(bestReferenceComplex.getName() + "\t");
			} else {
				sbGenes.append("NA\t");
			}
			final StringBuilder sbProteins = new StringBuilder();
			sbProteins.append(sbGenes.toString());
			// components genes
			sbGenes.append(predictedComplex.getComponentsGeneString("\t") + "\n");
			fwGenes.write(sbGenes.toString());
			// components proteins
			sbProteins.append(predictedComplex.getComponentsProteinsAccString("\t") + "\n");
			fwProteins.write(sbProteins.toString());

		}
		fwGenes.close();
		fwProteins.close();
		System.out.println("Predicted complexes written in file at: " + fileGenes.getAbsolutePath());
		System.out.println("Predicted complexes written in file at: " + fileProteins.getAbsolutePath());
	}

	private List<ProteinComplex> clusterInteractions(ClassificationResult result, double pValueCutOff)
			throws IOException {
		double clusterOneMinDensity = 0.4;
		if (optimizeClusterOneMinDensity) {
			clusterOneMinDensity = 0.0;
		}
		final double clusterOnePenalty = 2.9;
		double maxMatchingRatio = 0;
		List<ProteinComplex> finalClusters = null;
		int numIterationsWithNoImprovement = 0;
		while (clusterOneMinDensity < 1.0) {
			List<ProteinComplex> clusters = null;
			try {
				// OPTIMIZE this parameter to get maximum matching ratio between
				// the
				// known complexes and the resulting ones
				clusters = clusterInComplexes(result, pValueCutOff, ClassLabel.INTRA_COMPLEX, clusterOneMinDensity,
						clusterOnePenalty);
			} catch (final IllegalArgumentException e) {
				e.printStackTrace();
				log.warn("No clusters with pValueCutoff = " + pValueCutOff);
			}
			if (optimizeClusterOneMinDensity && clusters != null && !clusters.isEmpty()) {

				final double matchingRatio = ClusterEvaluation.getMaximumMatchingRatio(clusters,
						getReferenceDBSimplified().getProteinComplexes());
				log.info("MMR = " + matchingRatio + " with clusterOneMinDenbsity=" + clusterOneMinDensity);
				if (matchingRatio > maxMatchingRatio) {
					finalClusters = clusters;
					maxMatchingRatio = matchingRatio;
					log.info("Optimal MMR = " + matchingRatio + " found with clusterOneMinDenbsity="
							+ clusterOneMinDensity);
				} else {
					numIterationsWithNoImprovement++;
				}
			} else {
				numIterationsWithNoImprovement++;
				finalClusters = clusters;
			}
			if (!optimizeClusterOneMinDensity
					|| numIterationsWithNoImprovement == MAX_NUM_ITERATIONS_WITH_NO_IMPROVEMENT) {
				break;
			}
			clusterOneMinDensity = increase(clusterOneMinDensity);
		}

		if (finalClusters != null) {
			final double compositeScore = ClusterEvaluation.calculateCompositeScore(maxMatchingRatio, finalClusters,
					getReferenceDBSimplified().getProteinComplexes());
			log.info("COMPOSITE SCORE = " + compositeScore);
		}
		return finalClusters;
	}

	private void createCurvesImages(ClassificationResult result) throws IOException {
		final List<ProteinPairInteraction> interactions = result.getInteractions(0.0, ClassLabel.INTRA_COMPLEX);
		final double[][] prCurve = ClassificationEvaluation.getPRCurve(interactions);
		final File imageFile = new File(
				matrixFolder + File.separator + "PRCurve - " + experiment.getProjectName() + ".png");
		ImageGenerator.generateImage(imageFile, prCurve, experiment.getProjectName(), "Recall", "Precision");

		final double[][] rocCurve = ClassificationEvaluation.getROCCurve(interactions);
		final File imageFile2 = new File(
				matrixFolder + File.separator + "ROCCurve - " + experiment.getProjectName() + ".png");
		ImageGenerator.generateImage(imageFile2, rocCurve, experiment.getProjectName(), "False Positive Rate",
				"True Positive Rate");
	}

	private double increase(double clusterOneMinDensity) {
		return clusterOneMinDensity + 0.1;
	}

	private List<ProteinComplex> clusterInComplexes(ClassificationResult result, double pValueCutoff,
			ClassLabel classLabel, double clusterOneMinDensity, double clusterOnePenalty) {
		return new ClusterOneInterface().runClusterOne(result, pValueCutoff, classLabel, clusterOneMinDensity,
				clusterOnePenalty);
	}

	private double estimatePValueThreshold(ClassificationResult result, double precisionCutOff) {
		// sort by the probability of being intra-complex and calculate the
		// TP/(TP+FP)
		double pValueCutOff = 0.5;
		while (pValueCutOff <= 1.0) {
			final List<ProteinPairInteraction> interactions = result.getInteractions(pValueCutOff,
					ClassLabel.INTRA_COMPLEX);
			final double precision = ClassificationEvaluation.calculatePrecision(interactions);
			final double recall = ClassificationEvaluation.calculateRecall(interactions);
			final double f2Score = ClassificationEvaluation.calculateFMeasure(interactions);
			if (Double.isNaN(precision)) {
				pValueCutOff -= 0.02;
				break;
			}
			log.info(" pValue>=" + pValueCutOff + " -> precision=" + precision + ", recall=" + recall + ", f2score="
					+ f2Score);
			if (precision > precisionCutOff) {
				break;
			}
			pValueCutOff += 0.02;
			log.info("Now trying with p-Value < " + pValueCutOff);
		}
		pValueCutOff = Math.min(1.0, pValueCutOff);
		log.info("pValue cuttoff for getting a desired precision of " + precisionCutOff + " is = " + pValueCutOff);
		return pValueCutOff;
	}

	private void printOutSummary(PrintStream out, int numFP, int numNew, int numTP, double precisionCutOff) {
		final DecimalFormat formatter = new DecimalFormat("#.#");
		final int total = numFP + numNew + numTP;
		out.append("After applying threshold of precision <= " + formatter.format(precisionCutOff) + "\n");
		out.append("Number of known interactors: " + numTP + " (" + formatter.format(numTP * 100.0 / total) + ")\n");
		out.append("Number of new interactors: " + numNew + " (" + formatter.format(numNew * 100.0 / total) + ")\n");
		out.append("Number of false positive interactors: " + numFP + " (" + formatter.format(numFP * 100.0 / total)
				+ ")\n");
	}

	private ClassificationResult naiveBayesTraining(CoFractionationDataset dataset) throws Exception {

		if (machineLearningLibraryToUse == MLLIBRARY.SMILE) {
			throw new IllegalArgumentException("Matching learning library SMILE is not supported anymore");
			// final NaiveBayesSmile naiveBayesSmile = new NaiveBayesSmile(this,
			// matrixesByDistanceMapTraining,
			// proteinListTraining);
			// final MyClassifier classifier =
			// naiveBayesSmile.naiveBayesTraining();
			// final Instances instances = getInstances(matrixesByDistanceMap,
			// proteinList);
			// final ClassificationResult result =
			// classifier.evaluateClassifier(this, instances, proteinList);
			// return result;
		} else {
			// training

			final NaiveBayesWeka naiveBayesWeka = new NaiveBayesWeka(kFoldCrossValidation, minCorrelationForTraining,
					dataset.getTrainingDataset());
			final MyClassifier classifier = naiveBayesWeka.naiveBayesTraining();
			// make the classification in the whole dataset
			final ClassificationResult result = classifier.evaluateClassifier(dataset.getDataset());
			// // all
			// final File arffFile = new File(matrixFolder.getAbsolutePath() +
			// File.separator + projectName + "_all.arff");
			// NaiveBayesWeka.createArffFile(arffFile, matrixesByDistanceMap,
			// proteinList, -1.0, this, false);
			// final ArffLoader loader = new ArffLoader();
			// loader.setSource(arffFile);
			//
			// final Instances instances = loader.getDataSet();
			// final ClassificationResult result =
			// classifier.evaluateClassifier(this, instances, proteinList);
			return result;
		}

	}

	@Override
	public ClassLabel getClassLabel(Collection<String> proteins, Collection<String> proteins2) throws IOException {
		final ProteinComplexDB ref = getReferenceDBSimplified();
		final Set<ProteinComplex> proteinComplexes1 = ref.getProteinComplexesByProteins(proteins);
		final Set<ProteinComplex> proteinComplexes2 = ref.getProteinComplexesByProteins(proteins2);
		if (proteinComplexes1.isEmpty() || proteinComplexes2.isEmpty()) {
			// return ClassLabel.NOVEL_INTERACTORS;
			return null;
		}
		for (final ProteinComplex proteinComplex : proteinComplexes1) {
			if (proteinComplexes2.contains(proteinComplex)) {
				return ClassLabel.INTRA_COMPLEX;
			}
		}
		return ClassLabel.INTER_COMPLEX;
	}

	@Override
	public ClassLabel getClassLabel(String protein, String protein2) throws IOException {
		final ProteinComplexDB ref = getReferenceDBSimplified();
		final Set<ProteinComplex> proteinComplexes1 = ref.getProteinComplexesByProtein(protein);
		final Set<ProteinComplex> proteinComplexes2 = ref.getProteinComplexesByProtein(protein2);
		if (proteinComplexes1.isEmpty() || proteinComplexes2.isEmpty()) {
			// return ClassLabel.NOVEL_INTERACTORS;
			return null;
		}
		for (final ProteinComplex proteinComplex : proteinComplexes1) {
			if (proteinComplexes2.contains(proteinComplex)) {
				return ClassLabel.INTRA_COMPLEX;
			}
		}
		return ClassLabel.INTER_COMPLEX;
	}

	private void loadDistancesMatrixes(Map<DistanceMeasure, NLMatrix> distanceMap) throws IOException {
		distanceMap.clear();
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			final File matrixFile = getOutputFileNameForMatrix(distance);
			if (matrixFile.exists()) {
				try {
					final NLMatrix matrix = loadDistanceMatrix(matrixFile);
					distanceMap.put(distance, matrix);
				} catch (final IllegalArgumentException e) {
					log.warn(e.getMessage());
					log.warn("Skipping loading of distance " + distance.name() + " from matrix in file "
							+ matrixFile.getAbsolutePath());
				}
			}
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
				for (int col = row + 1; col < split.length; col++) {
					final double num = Double.valueOf(split[col]);
					if (Double.isNaN(num)) {
						throw new IllegalArgumentException(
								"There is a NaN value at matrix " + matrixFile.getAbsolutePath());
					}
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
		// in case of pvalues, make the log first
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			if (distance == DistanceMeasure.CORRELATION_PVALUE) {
				makelog2Matrix(distanceMap.get(distance));
			}
		}
		// standardize values in matrixes
		// making mean=0 and standard deviation 1
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			standardizeMatrix(distanceMap.get(distance));
		}
		log.info(DistanceMeasure.values().length + " matrixes pre-processed");
	}

	private void calculateDistanceMeasurements(List<String> proteinList, Map<DistanceMeasure, NLMatrix> distanceMap,
			DataType dataType, List<Boolean> calculateDistance) throws IOException, InvalidConfigurationException {
		final FileWriter logFileWriter = new FileWriter(logFile, false);
		final long total = proteinList.size() * (proteinList.size() - 1l) / 2;
		final ProgressCounter counter = new ProgressCounter(total, ProgressPrintingType.PERCENTAGE_STEPS, 0);
		counter.setShowRemainingTime(true);
		for (int i = 0; i < proteinList.size(); i++) {
			final String protein1 = proteinList.get(i);
			for (int j = i + 1; j < proteinList.size(); j++) {
				final String protein2 = proteinList.get(j);
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info(printIfNecessary);
					if (!fittingTimes.isEmpty()) {
						final double averageTime = Maths.mean(fittingTimes);
						log.info(DatesUtil.getDescriptiveTimeFromMillisecs(averageTime) + " in average per protein");
					}
				}
				for (final DistanceMeasure distance : DistanceMeasure.values()) {
					if (calculateDistance.get(distance.ordinal())) {
						switch (distance) {
						case PEARSON_CORRELATION_COEFFICIENT:
							final double correlationCoefficient = PearsonCorrelation
									.getPearsonCorrelationCoefficient(experiment, protein1, protein2, dataType);
							distanceMap.get(distance).set(i, j, correlationCoefficient);
							log.debug("correlationCoefficient" + "\t" + correlationCoefficient);
							break;
						case CORRELATION_PVALUE:
							final double pvalue = PearsonCorrelation.getPearsonCorrelationPValue(experiment, protein1,
									protein2, dataType);
							distanceMap.get(distance).set(i, j, pvalue);
							log.debug("pvalue" + "\t" + pvalue);
							break;
						case EUCLIDEAN_DISTANCE:
							final double euDistance = EuclideanDistanceCalculator.getEuclideanDistance(experiment,
									protein1, protein2, dataType);
							distanceMap.get(distance).set(i, j, euDistance);
							log.debug("euDistance" + "\t" + euDistance);
							break;
						case PEAK_LOCATION:
							final double peakLocation = getPeakLocation(experiment, protein1, protein2, dataType);
							log.debug("peakLocation" + "\t" + peakLocation);
							distanceMap.get(distance).set(i, j, peakLocation);
							break;
						case COAPEX_SCORE:
							final double coApexScore = getCoApexScore(experiment, protein1, protein2, dataType,
									logFileWriter, logFittings);
							distanceMap.get(distance).set(i, j, coApexScore);
							log.debug("coApexScore" + "\t" + coApexScore);
							break;
						case APEX_SCORE:
							final double apexScore = ApexScore.getApexScore(experiment, protein1, protein2, dataType);
							distanceMap.get(distance).set(i, j, apexScore);
							log.debug("coApexScore" + "\t" + apexScore);
							break;
						// case BAYES_CORRELATION:
						// final double bayesCorrelation =
						// getBayesCorrelation(experiment, protein1, protein2,
						// dataType);
						// distanceMap.get(distance).set(i, j,
						// bayesCorrelation);
						// break;
						case JACCARD_SCORE:
							final double jaccardScore = JaccardScore.getJaccardScore(experiment, protein1, protein2);
							distanceMap.get(distance).set(i, j, jaccardScore);
							break;
						case MUTUAL_INFORMATION:
							final double mutualInformation = MutualInformation.getMutualInformation(experiment,
									protein1, protein2);
							distanceMap.get(distance).set(i, j, mutualInformation);
							break;
						case PEARSON_CORRELATION_COEFFICIENT_PLUS_NOISE:
							final double pccn = PearsonCorrelationPlusNoise.getInstance(experiment)
									.getPearsonCorrelationPlusNoise(protein1, protein2);
							distanceMap.get(distance).set(i, j, pccn);
							break;
						case WEIGTHED_CROSS_CORRELATION:
							final double wcc = WeightedCrossCorrelation.getWeightedCrossCorrelation(experiment,
									protein1, protein2, dataType);
							distanceMap.get(distance).set(i, j, wcc);
							break;
						default:
							break;
						}
					}
				}

			}
		}
	}

	private File getOutputFileNameForMatrix(DistanceMeasure distanceMeasure) {
		return new File(matrixFolder.getAbsolutePath() + File.separator + projectName + "_" + distanceMeasure.name()
				+ ".matrix");
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
			referenceDB = ProteinComplexAnalyzer.getDBs().get(0);

		}
		return referenceDB;

	}

	private ProteinComplexDB getReferenceDBSimplified() throws IOException {

		if (referenceDBSimplified == null) {
			referenceDBSimplified = ProteinComplexAnalyzer.getDBs(true, maxOverlapScoreInReferenceSet,
					minComplexSizeInReferenceSetForLearning, maxComplexSizeInReferenceSetForLearning).get(0);

		}
		return referenceDBSimplified;

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
		final TDoubleList profile1 = exp.getNormalizedElutionProfile(acc1, dataType);
		final TDoubleList profile2 = exp.getNormalizedElutionProfile(acc2, dataType);

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

	private void printFitting(List<MyGaussianFit> fits, TDoubleList profile, String acc, FileWriter fw)
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

	private List<MyGaussianFit> fitToGaussiansWithGeneticAlgorithm(TDoubleList profile, String acc, TIntList spcList)
			throws InvalidConfigurationException {
		if (fitsByProteins.containsKey(acc)) {
			return fitsByProteins.get(acc);
		}
		final long t0 = System.currentTimeMillis();
		// fill missing values with the average of neighbor points
		TDoubleList profileWithImputedMissingValues = null;
		if (averageMissingValues) {
			profileWithImputedMissingValues = imputeMissingValuesOfProfile(profile);
		} else {
			profileWithImputedMissingValues = profile;
		}
		// smooth profile by a sliding average of width 3
		final TDoubleList profileSmoothed = smoothProfile(profileWithImputedMissingValues, PROFILE_SMOOTH_WIDTH);

		log.info("Fitting profile for protein " + acc);

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
		window.setAverageFittingTime(Maths.mean(fittingTimes));
		return ret;
	}

	private ModelPlot launchFitting(String acc, int numGaussians, TIntList spcProfile, TDoubleList rawProfile,
			TDoubleList processedProfile) throws InvalidConfigurationException {
		final MultipleGaussianFitModel model = new MultipleGaussianFitModel(acc, numGaussians);

		model.setProcessedProfile(processedProfile);
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
	private TDoubleList imputeMissingValuesOfProfile(TDoubleList profile) {
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

		final double distance = EuclideanDistanceCalculator.compute(array1.toArray(), array2.toArray());
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
	private TIntList getFractionsWithMaximums(TDoubleList profile1) {
		final TIntList ret = new TIntArrayList();
		final double max = profile1.max();

		for (int fraction = 0; fraction < profile1.size(); fraction++) {
			if (Double.compare(profile1.get(fraction), max) == 0) {
				ret.add(fraction);
			}
		}
		return ret;
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
		List<MyGaussianFit> fits1 = null;
		final FittingCache fittingCache = FittingCache.getInstance(exp);
		if (fittingCache.containsFittingsForProtein(acc1)) {
			fits1 = fittingCache.getFittingsForProtein(acc1);
		} else {
			final TDoubleList profile1 = exp.getNormalizedElutionProfile(acc1, dataType);
			final TIntList spcprofile1 = getDataProfileIntegers(exp, acc1, DataType.SPC);
			fits1 = fitToGaussiansWithGeneticAlgorithm(profile1, acc1, spcprofile1);
			fittingCache.put(acc1, fits1);
			if (logFittings) {
				printFitting(fits1, profile1, acc1, fw);
			}
		}
		if (fits1.isEmpty()) {
			throw new IllegalArgumentException("This shoudn't happen");
		}
		List<MyGaussianFit> fits2 = null;
		if (fittingCache.containsFittingsForProtein(acc2)) {
			fits2 = fittingCache.getFittingsForProtein(acc2);
		} else {
			final TDoubleList profile2 = exp.getNormalizedElutionProfile(acc2, dataType);
			final TIntList spcprofile2 = getDataProfileIntegers(exp, acc2, DataType.SPC);
			fits2 = fitToGaussiansWithGeneticAlgorithm(profile2, acc2, spcprofile2);
			fittingCache.put(acc2, fits2);
			if (logFittings) {
				printFitting(fits2, profile2, acc2, fw);
			}
		}
		if (fits2.isEmpty()) {
			throw new IllegalArgumentException("This shoudn't happen");
		}
		final double minimumDistance = getDistanceBetweenClosestGaussians(fits1, fits2);
		return minimumDistance;
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

	/**
	 * Take all non Nan values from the matrix and normalize them so that the
	 * mean is 0 and the standard deviation is 1
	 * 
	 * @param nlMatrix
	 */
	public void makelog2Matrix(NLMatrix matrix) {
		final long t1 = System.currentTimeMillis();
		log.info("Making log2 from  matrix of " + matrix.nrows() + "x" + matrix.ncols());
		final int cols = matrix.ncols();
		final int rows = matrix.nrows();

		final TDoubleArrayList nums = new TDoubleArrayList();
		final TDoubleArrayList numsLog = new TDoubleArrayList();
		for (int row = 0; row < rows; row++) {
			for (int col = 0; col < cols; col++) {
				double num = matrix.get(row, col);
				if (!Double.isNaN(num)) {
					nums.add(num);
					if (num == 0.0) {
						num = Double.MIN_VALUE;
					}
					final double log2 = Maths.log(num, 2);
					numsLog.add(log2);
					matrix.set(row, col, log2);
				}
			}
		}
		final long t2 = System.currentTimeMillis();
		log.info("Matrix of " + matrix.nrows() + "x" + matrix.ncols() + " are not log 2 values. done in "
				+ DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));
		log.info("Before: min=" + nums.min() + " max=" + nums.max() + " and after log2: min=" + numsLog.min() + " max="
				+ numsLog.max());
	}
}
