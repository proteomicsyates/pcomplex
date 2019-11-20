package edu.scripps.yates.pcomplex.distances;

import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import smile.netlib.NLMatrix;
import smile.stat.distribution.PoissonDistribution;

public class PearsonCorrelationPlusNoise {
	private final static Logger log = Logger.getLogger(PearsonCorrelationPlusNoise.class);
	private static int N = 5;
	private final NLMatrix matrix;
	private final SeparationExperiment exp;
	private static final Map<SeparationExperiment, PearsonCorrelationPlusNoise> pearsonCorrelations = new THashMap<SeparationExperiment, PearsonCorrelationPlusNoise>();

	private PearsonCorrelationPlusNoise(SeparationExperiment exp) {
		this.exp = exp;
		matrix = calculatePearsonCorrelationPlusNoiseMatrix(exp);
	}

	public static PearsonCorrelationPlusNoise getInstance(SeparationExperiment exp) {
		if (!pearsonCorrelations.containsKey(exp)) {
			pearsonCorrelations.put(exp, new PearsonCorrelationPlusNoise(exp));
		}
		return pearsonCorrelations.get(exp);
	}

	public double getPearsonCorrelationPlusNoise(String protein1, String protein2) {
		final int i = exp.getIndexForProtein(protein1);
		final int j = exp.getIndexForProtein(protein2);
		return matrix.get(i, j);
	}

	private static NLMatrix calculatePearsonCorrelationPlusNoiseMatrix(SeparationExperiment exp) {
		log.info("Calculating matrix for pearson correlation with modeled noise");
		final ProgressCounter counter = new ProgressCounter(N, ProgressPrintingType.PERCENTAGE_STEPS, 0);
		counter.setShowRemainingTime(true);
		final NLMatrix ret = new NLMatrix(exp.getNumProteins(), exp.getNumProteins());
		for (int step = 1; step <= N; step++) {
			final NLMatrix profileMatrixWithNoise = createProfileMatrixWithNoise(exp);
			final double[][] array = profileMatrixWithNoise.array();
			// that matrix has as many rows as proteins in the experiment and as
			// many columns as fractions.
			final int nrows = profileMatrixWithNoise.nrows();
			for (int i = 0; i < nrows; i++) {
				// protein i

				final double[] profile1 = array[i];
				for (int j = i + 1; j < nrows; j++) {
					final double[] profile2 = array[j];
					// calculate the pearson correlation of both proteins
					final double pcoef = PearsonCorrelation.correlationCoefficient(profile1, profile2);
					if (step == 1) {
						ret.set(i, j, pcoef);
					} else {
						final double previousPCoeff = ret.get(i, j);
						final double mean = ((previousPCoeff * (step - 1)) + pcoef) / step;
						ret.set(i, j, mean);
					}
				}
			}
			counter.increment();
			final String progress = counter.printIfNecessary();
			if (!"".equals(progress)) {
				log.info(progress);
			}
		}
		return ret;
	}

	public static NLMatrix createProfileMatrixWithNoise(SeparationExperiment exp) {
		log.info("Creating matrix with noise");
		final double noiseConstant = 1.0 / exp.getFractions().size();
		final NLMatrix profileWithNoise = new NLMatrix(exp.getNumProteins(), exp.getFractions().size(), noiseConstant);

		final TIntList spcs = new TIntArrayList();
		for (int i = 0; i < exp.getNumProteins(); i++) {
			final String acc = exp.getProteinACCList().get(i);
			for (int j = 0; j < exp.getFractions().size(); j++) {
				final Fraction fraction = exp.getFractions().get(j);
				final Protein protein = fraction.getProteinByAcc(acc);
				if (protein != null) {
					spcs.add(Double.valueOf(protein.getData(DataType.SPC)).intValue());
				} else {
					spcs.add(0);
				}
			}
		}
		// use the poisson function
		final PoissonDistribution poisson = new PoissonDistribution(spcs.toArray());
		log.info("Poisson distribution estimated with lambda=" + poisson.getLambda());
		for (int i = 0; i < exp.getNumProteins(); i++) {
			final String acc = exp.getProteinACCList().get(i);
			for (int j = 0; j < exp.getFractions().size(); j++) {
				final Fraction fraction = exp.getFractions().get(j);
				final Protein protein = fraction.getProteinByAcc(acc);
				final double rand = poisson.rand();
				if (protein != null) {
					final double spc = protein.getData(DataType.SPC);

					final double c = spc + rand;
					profileWithNoise.set(i, j, c);
				} else {
					final double c = rand;
					profileWithNoise.set(i, j, c);
				}
			}
		}
		log.info("Normalizing matrix");
		// now normalize by rows
		for (int i = 0; i < exp.getNumProteins(); i++) {
			final TDoubleArrayList nums = new TDoubleArrayList();
			for (int j = 0; j < exp.getFractions().size(); j++) {
				nums.add(profileWithNoise.get(i, j));
			}
			for (int j = 0; j < exp.getFractions().size(); j++) {
				profileWithNoise.set(i, j, profileWithNoise.get(i, j) / nums.sum());
			}
		}
		// now normalize by columns
		for (int j = 0; j < exp.getFractions().size(); j++) {
			final TDoubleArrayList nums = new TDoubleArrayList();
			for (int i = 0; i < exp.getNumProteins(); i++) {
				nums.add(profileWithNoise.get(i, j));
			}
			for (int i = 0; i < exp.getNumProteins(); i++) {
				profileWithNoise.set(i, j, profileWithNoise.get(i, j) / nums.sum());
			}
		}
		log.info("Matrix with noise ready");
		return profileWithNoise;
	}
}
