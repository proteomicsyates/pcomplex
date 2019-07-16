package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.Logger;
import org.jgap.Configuration;
import org.jgap.Gene;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.jgap.NaturalSelector;
import org.jgap.Population;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;

public class GaussianNaturalSelector extends NaturalSelector {
	private final static Logger log = Logger.getLogger(GaussianNaturalSelector.class);
	/**
	 * 
	 */
	private static final long serialVersionUID = 4101939976816364683L;

	public GaussianNaturalSelector(Configuration conf) throws InvalidConfigurationException {
		super(conf);
	}

	private final List<IChromosome> chromoromes = new ArrayList<IChromosome>();

	@Override
	public void empty() {
		chromoromes.clear();
	}

	@Override
	public boolean returnsUniqueChromosomes() {

		return false;
	}

	@Override
	protected void add(IChromosome a_chromosomeToAdd) {
		chromoromes.add(a_chromosomeToAdd);

	}

	@Override
	public void select(int a_howManyToSelect, Population a_from_population, Population a_to_population) {

		final List<IChromosome> chromosomes = a_from_population.getChromosomes();
		// sort by fitness
		Collections.sort(chromosomes, new Comparator<IChromosome>() {

			@Override
			public int compare(IChromosome o1, IChromosome o2) {
				return Double.compare(o2.getFitnessValue(), o1.getFitnessValue());
			}
		});
		boolean minimumSelected = false;
		;
		int numSelected = 0;
		for (final IChromosome iChromosome : chromosomes) {

			if (numSelected == a_howManyToSelect) {
				minimumSelected = true;
				break;
			}
			log.debug(iChromosome.getFitnessValue());
			final boolean valid = isValid(iChromosome);
			if (valid) {
				numSelected++;
			}
			if (valid) {
				a_to_population.addChromosome(iChromosome);
			}
		}

	}

	private boolean isValid(IChromosome iChromosome) {
		final Gene[] genes = iChromosome.getGenes();
		final List<MyGaussianFit> gaussians = FittingUtil.getGaussiansFromGenes(genes);
		if (thereAreSubsetGaussians(gaussians)) {
			return false;
		}
		return true;
	}

	private boolean thereAreSubsetGaussians(List<MyGaussianFit> gaussians) {
		for (int i = 0; i < gaussians.size(); i++) {
			final MyGaussianFit gaussian1 = gaussians.get(i);
			for (int j = i + 1; j < gaussians.size(); j++) {
				final MyGaussianFit gaussian2 = gaussians.get(j);

				final boolean isSubset = isSubsetGaussians(gaussian1, gaussian2);
				if (isSubset) {
					// System.out.print("Gaussian: ");
					// FittingUtil.printGaussian(System.out, gaussian1);
					// System.out.print(" is subset of gaussian:");
					// FittingUtil.printGaussian(System.out, gaussian2);
					return true;
				}
				final boolean isSubset2 = isSubsetGaussians(gaussian2, gaussian1);
				if (isSubset2) {
					// System.out.print("Gaussian: ");
					// FittingUtil.printGaussian(System.out, gaussian2);
					// System.out.print(" is subset of gaussian:");
					// FittingUtil.printGaussian(System.out, gaussian1);
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Returns true if gaussian1 is subset of gaussian2
	 * 
	 * @param gaussian1
	 * @param gaussian2
	 * @return
	 */
	private boolean isSubsetGaussians(MyGaussianFit gaussian1, MyGaussianFit gaussian2) {
		final double mean_1 = FittingUtil.getMean(gaussian1.getFittedParameters());
		final double mean_2 = FittingUtil.getMean(gaussian2.getFittedParameters());

		if (mean_1 < mean_2) {
			final double minus2Sigma_2 = FittingUtil.getXSigma(gaussian2.getFittedParameters(), -2);
			if (minus2Sigma_2 < mean_1) {
				return true;
			}
		} else {
			final double plus2Sigma_2 = FittingUtil.getXSigma(gaussian2.getFittedParameters(), 2);
			if (plus2Sigma_2 > mean_1) {
				return true;
			}
		}

		return false;
	}

}
