package edu.scripps.yates.pcomplex.model;

import java.util.Set;

public class ProteinComplexExistenceCriteria {
	private final static double percentageInSmallComplexes = 50; // 50%
	private final static double percentageInBigComplexes = 30;// 30%

	public boolean considersExisting(ProteinComplex proteinComplex, Set<String> identifiedProteins) {
		if (proteinComplex.getComponents().size() < 5) {
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

}
