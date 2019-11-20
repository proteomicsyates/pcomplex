package edu.scripps.yates.pcomplex.model;

import java.util.Set;

public class ProteinComplexExistenceCriteria {
	public static final double defaultPercentageInSmallComplexes = 50; // 50%
	public static final double defaultPercentageInBigComplexes = 30;// 30%
	public static final int default_MAX_COMPONENTS_TO_BE_SMALL = 5;
	private double percentageInSmallComplexes = defaultPercentageInSmallComplexes; // 50%
	private double percentageInBigComplexes = defaultPercentageInBigComplexes;// 30%
	private int maxComponentsToBeSmall = default_MAX_COMPONENTS_TO_BE_SMALL;

	public boolean considersExisting(ProteinComplex proteinComplex, Set<String> identifiedProteins) {
		if (proteinComplex.getComponents().size() < maxComponentsToBeSmall) {
			return considersExistingInSmallProteinComplexesCriteria(proteinComplex, identifiedProteins);
		} else {
			return considersExistingInBigProteinComplexesCriteria(proteinComplex, identifiedProteins);
		}
	}

	private boolean considersExistingInBigProteinComplexesCriteria(ProteinComplex proteinComplex,
			Set<String> proteinsAndGenes) {
		final double percentageOfUnits = getPercentageOfUnits(proteinComplex, proteinsAndGenes);
		if (percentageOfUnits >= percentageInBigComplexes) {
			return true;
		}
		return false;
	}

	private double getPercentageOfUnits(ProteinComplex proteinComplex, Set<String> proteinsAndGenes) {
		final int totalUnitsInComplex = proteinComplex.getComponents().size();
		final int unitsCovered = proteinComplex.getComponentsInProteinsAndGenes(proteinsAndGenes).size();
		if (unitsCovered == 0) {
			return 0.0;
		}
		return unitsCovered * 100.0 / totalUnitsInComplex;
	}

	/**
	 * 
	 * @param proteinComplex
	 * @param identifiedProteinsKeys
	 * @return
	 */
	private boolean considersExistingInSmallProteinComplexesCriteria(ProteinComplex proteinComplex,
			Set<String> proteinsAndGenes) {
		final double percentageOfUnits = getPercentageOfUnits(proteinComplex, proteinsAndGenes);
		if (percentageOfUnits >= percentageInSmallComplexes) {
			return true;
		}
		return false;
	}

	public double getPercentageInSmallComplexes() {
		return percentageInSmallComplexes;
	}

	public void setPercentageInSmallComplexes(double percentageInSmallComplexes) {
		this.percentageInSmallComplexes = percentageInSmallComplexes;
	}

	public double getPercentageInBigComplexes() {
		return percentageInBigComplexes;
	}

	public void setPercentageInBigComplexes(double percentageInBigComplexes) {
		this.percentageInBigComplexes = percentageInBigComplexes;
	}

	public int getMaxComponentsToBeSmall() {
		return maxComponentsToBeSmall;
	}

	public void setMaxComponentsToBeSmall(int maxComponentsToBeSmall) {
		this.maxComponentsToBeSmall = maxComponentsToBeSmall;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append("A protein complex with at most " + this.maxComponentsToBeSmall
				+ " components will be present if at least " + this.percentageInSmallComplexes
				+ "% of their components have been detected.\n");
		sb.append("A protein complex with " + (this.maxComponentsToBeSmall + 1)
				+ " or more components will be present if at least " + this.percentageInBigComplexes
				+ "% of their components have been detected.\n");
		return sb.toString();
	}
}
