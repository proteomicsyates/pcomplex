package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.util.List;

import org.apache.log4j.Logger;
import org.jgap.Chromosome;
import org.jgap.Configuration;
import org.jgap.Gene;
import org.jgap.Genotype;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.jgap.NaturalSelector;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.DoubleGene;

import edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian.AbstractMultipleGaussianFit;
import edu.scripps.yates.pcomplex.cofractionation.fitting.gaussian.MultipleGaussianFit;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;

/**
 * tries to fit using a genetic alg. the ratio distribution to a mixture of
 * gaussians
 *
 * @author diego and Salva
 *
 */
public class MultipleGaussianFitModel extends AbstractMultipleFitModel {
	private final static Logger log = Logger.getLogger(MultipleGaussianFitModel.class);

	private static final int DEFAULT_MAX_ITERATIONS = 400;
	private static final int DEFAULT_POPULATION_SIZE = 500;

	public static final String TASK_PROGRESS = "progress";
	public static final String FITTING_RESULT = "fitting result";
	public static final String FITTING_RESULT_OPTIMAL = "fittin result optimal";

	private double fitnessValue;

	private double rSquare;

	public MultipleGaussianFitModel(String acc, int numGaussians) {
		super(acc, null, numGaussians);
	}

	@Override
	public ModelPlot runFit() throws InvalidConfigurationException {
		Configuration.reset();
		conf = new DefaultConfiguration();
		conf.setFitnessFunction(optimizationFunction);
		conf.setPopulationSize(DEFAULT_POPULATION_SIZE);
		log.debug("Starting fit with " + numGaussians + " gaussians");

		final double maxy = getMaxY();
		final List<double[]> experimentalData = getExperimentalData();
		final double maxx = experimentalData.get(experimentalData.size() - 1)[0];
		final double minx = experimentalData.get(0)[0];

		/*
		 * percentMaxY and the removal of maxx calc were added on 6-22-17 The
		 * second gaussian was being dropped to a height of 0, so we decided to
		 * fix the minimum height of the second gaussian to be at a minimum % of
		 * the max. Not sure what was wrong with maxx, but the second gaussian
		 * would have a mean where no datapoints lie.
		 */
		final double percentMaxY = 0.05;
		final double miny = (percentMaxY * maxy);
		log.debug("The minimum height of a gaussian would be " + miny + " wich is a " + (percentMaxY * 100)
				+ "% of the maximum Y value (" + maxy + ")");
		// maxx = Math.abs(Math.max(Math.abs(maxx), Math.abs(minx)));

		optimizationFunction.setExperimentalData(experimentalData);

		final DoubleGene[] genes = new DoubleGene[3 * numGaussians];
		for (int numGaussian = 0; numGaussian < numGaussians; numGaussian++) {

			if (fixedHeight[numGaussian] != null) {
				genes[numGaussian * 3] = new DoubleGene(conf, fixedHeight[numGaussian], fixedHeight[numGaussian]); // A
			} else {
				genes[numGaussian * 3] = new DoubleGene(conf, miny, maxy + 0.3 * maxy); // A
			}
			if (fixedAvg[numGaussian] != null) {
				genes[numGaussian * 3 + 1] = new DoubleGene(conf, fixedAvg[numGaussian], fixedAvg[numGaussian]); // avg
			} else {
				genes[numGaussian * 3 + 1] = new DoubleGene(conf, 0,
						// minx - 0.2 * maxx, maxx + 0.2 * maxx
						maxx + 1); // avg
			}
			if (fixedStd[numGaussian] != null) {
				genes[numGaussian * 3 + 2] = new DoubleGene(conf, fixedStd[numGaussian], fixedStd[numGaussian]); // std

			} else {
				genes[numGaussian * 3 + 2] = new DoubleGene(conf, 1, FittingUtil.MAX_GAUSSIAN_WIDTH / 2); // std
			}

		}
		final Chromosome sampleChromosome = new Chromosome(conf, genes);

		conf.setSampleChromosome(sampleChromosome);

		conf.setPopulationSize(1000);
		final NaturalSelector naturalSelector = new GaussianNaturalSelector(conf);
		conf.addNaturalSelector(naturalSelector, true);

		final Genotype genotype = Genotype.randomInitialGenotype(conf);

		progress = 0;
		final ProgressCounter counter = new ProgressCounter(getMaxIterations(), ProgressPrintingType.PERCENTAGE_STEPS,
				0);
		counter.setShowRemainingTime(true);
		double previousFitnessValue = 0;
		int roundsWithNoChange = 0;
		for (int i = 0; i < getMaxIterations(); i++) {
			genotype.evolve();
			counter.increment();
			progress++;
			final String percentage = counter.printIfNecessary();
			if (!"".equals(percentage)) {
				// log.info(percentage);
			}

			final IChromosome bestSolutionSoFar = genotype.getFittestChromosome();
			final Gene[] best = bestSolutionSoFar.getGenes();

			fittedGaussians = FittingUtil.getGaussiansFromGenes(best);

			// modelPlot = new ModelPlot(experimentalData, this);
			fitnessValue = bestSolutionSoFar.getFitnessValue();
			rSquare = FittingUtil.calculateR2(fittedGaussians, getExperimentalData());
			// modelPlot.setrSquared(new double[] { fitnessValue,
			// leftFlankRSquared });

			final double changeValue = Math.abs(fitnessValue - previousFitnessValue);
			log.debug("Iteration change=" + changeValue);
			if (changeValue < CHANGE_FITNESS_THRESHOLD) {
				roundsWithNoChange++;
				if (roundsWithNoChange == MAX_ROUNDS_WITH_NO_CHANGE) {
					log.debug("Algorithm break");
					break;
				}
			} else {
				roundsWithNoChange = 0;
			}
			previousFitnessValue = fitnessValue;

		}

		// swingWorker.firePropertyChange(StatCalculator.FITTING_RESULT, null,
		// new Pair<Long, ModelPlot>(taskID, modelPlot));
		// swingWorker.firePropertyChange(CoPIT2.TASK_NAME, null, "Fitting
		// done.");

		log.debug("Returning optimized values  ");

		final ModelPlot model = new ModelPlot(acc, transformExperimentalDataForModelling(getSPCProfile()),
				transformExperimentalDataForModelling(getRawProfile()), getExperimentalData(), this,
				genotype.getPopulation());
		return model;
	}

	public double getModelB(int gaussianIndex) {
		return conf.getSampleChromosome().getGene(gaussianIndex * 3).getEnergy();
	}

	public double getModelAvgSig(int gaussianIndex) {
		return conf.getSampleChromosome().getGene(gaussianIndex * 3 + 1).getEnergy();
	}

	public double getModelStdSig(int gaussianIndex) {
		return conf.getSampleChromosome().getGene(gaussianIndex * 3 + 2).getEnergy();
	}

	@Override
	public AbstractMultipleGaussianFit createOptimizationFunction(Double maxLimit) {
		return new MultipleGaussianFit(maxLimit);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * edu.scripps.copit.statistical.fitting.AbstractFitModel#getMaxIterations()
	 */
	@Override
	public int getMaxIterations() {
		return DEFAULT_MAX_ITERATIONS;
	}

	@Override
	public double getFitnessValue() {
		return fitnessValue;
	}

	@Override
	public double getRSquare() {
		return rSquare;
	}

}
