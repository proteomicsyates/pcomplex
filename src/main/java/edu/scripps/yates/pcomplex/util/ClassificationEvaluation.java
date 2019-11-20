package edu.scripps.yates.pcomplex.util;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;

public class ClassificationEvaluation {
	private final static Logger log = Logger.getLogger(ClassificationEvaluation.class);

	public static double calculatePrecision(List<ProteinPairInteraction> predictedInteractions) {
		int numTP = 0;
		int numFP = 0;

		for (final ProteinPairInteraction ppi : predictedInteractions) {
			if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
				if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
					numTP++;
				} else if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
					numFP++;
				}
			}
		}
		final double precision = 1.0 * numTP / (numTP + numFP);
		// printOutSummary(System.out, numFP, numNew, numTP, precisionCutOff);
		return precision;
	}

	public static double calculateRecall(List<ProteinPairInteraction> predictedInteractions) {
		int numTP = 0;
		int numFN = 0;

		for (final ProteinPairInteraction ppi : predictedInteractions) {
			if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
				if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
					numTP++;
				}
			} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
				// predicted as false
				if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
					numFN++;
				}
			}
		}
		final double recall = 1.0 * numTP / (numTP + numFN);
		// printOutSummary(System.out, numFP, numNew, numTP, precisionCutOff);
		return recall;
	}

	private static double calculateFalsePositiveRate(List<ProteinPairInteraction> interactions) {
		int numFP = 0;
		int numTN = 0;
		for (final ProteinPairInteraction ppi : interactions) {
			if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
				if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
					numFP++;
				}
			} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
				// predicted as false
				if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
					numTN++;
				}
			}
		}
		final double fpr = 1.0 * numFP / (numFP + numTN);
		return fpr;
	}

	public static double calculateFMeasure(List<ProteinPairInteraction> predictedInteractors) {
		final double precision = calculatePrecision(predictedInteractors);
		final double recall = calculateRecall(predictedInteractors);
		final double fMeasure = 2.0 * (precision * recall) / (precision + recall);
		return fMeasure;
	}

	public static double[][] getROCCurve(List<ProteinPairInteraction> interactors) {
		log.info("Creating ROC curve");
		Collections.sort(interactors, new Comparator<ProteinPairInteraction>() {

			@Override
			public int compare(ProteinPairInteraction o1, ProteinPairInteraction o2) {
				return Double.compare(o1.getProbability(ClassLabel.INTRA_COMPLEX),
						o2.getProbability(ClassLabel.INTRA_COMPLEX));
			}
		});
		final double[][] ret = new double[interactors.size()][2];

		int numTP = 0;
		int numFP = 0;
		int numTN = 0;
		int numFN = 0;
		for (final ProteinPairInteraction ppi : interactors) {
			if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
				if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
					numTP++;
				} else if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
					numFP++;
				}
			} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
				// predicted as false
				if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
					numTN++;
				}
				// predicted as false
				else if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
					numFN++;
				}
			}
		}

		double fpr = 1.0 * numFP / (numFP + numTN);
		double recall = 1.0 * numTP / (numTP + numFN);
		int i = 0;
		ret[i][0] = fpr;
		ret[i][1] = recall;
		// now the first are the best ones
		for (i = 1; i < interactors.size(); i++) {
			final ProteinPairInteraction ppi = interactors.get(i - 1);
			if (ppi.getTrueLabel() != null) {
				if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
					if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numTP--;
					} else if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numFP--;
					}
				} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
					// predicted as false
					if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numTN--;
					} else if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numFN--;
					}
				}
				fpr = 1.0 * numFP / (numFP + numTN);
				recall = 1.0 * numTP / (numTP + numFN);
				if (!Double.isNaN(fpr) && !Double.isNaN(recall)) {
					ret[i][0] = fpr;
					ret[i][1] = recall;
				}
			}

		}
		log.info("ROC curve created");
		return ret;

	}

	public static double[][] getPRCurve(List<ProteinPairInteraction> interactors) {
		log.info("Creating PR curve ");
		Collections.sort(interactors, new Comparator<ProteinPairInteraction>() {

			@Override
			public int compare(ProteinPairInteraction o1, ProteinPairInteraction o2) {
				return Double.compare(o1.getProbability(ClassLabel.INTRA_COMPLEX),
						o2.getProbability(ClassLabel.INTRA_COMPLEX));
			}
		});
		// the first ones are the ones with INTRA_COMPLEX high probabilities

		final double[][] ret = new double[interactors.size()][2];

		int numTP = 0;
		int numFP = 0;
		int numTN = 0;
		int numFN = 0;
		for (final ProteinPairInteraction ppi : interactors) {
			if (ppi.getTrueLabel() != null) {
				if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
					if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numTP++;
					} else if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numFP++;
					}
				} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
					// predicted as false
					if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numTN++;
					}
					// predicted as false
					else if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numFN++;
					}
				}
			}
		}

		double precision = 1.0 * numTP / (numTP + numFP);
		double recall = 1.0 * numTP / (numTP + numFN);
		int i = 0;
		ret[i][0] = recall;
		ret[i][1] = precision;
		// now the first are the best ones
		for (i = 1; i < interactors.size(); i++) {
			final ProteinPairInteraction ppi = interactors.get(i - 1);
			if (ppi.getTrueLabel() != null) {
				if (ppi.getPredictedLabel() == ClassLabel.INTRA_COMPLEX) {
					if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numTP--;
					} else if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numFP--;
					}
				} else if (ppi.getPredictedLabel() == ClassLabel.INTER_COMPLEX) {
					// predicted as false
					if (ppi.getTrueLabel() == ClassLabel.INTER_COMPLEX) {
						numTN--;
					} else if (ppi.getTrueLabel() == ClassLabel.INTRA_COMPLEX) {
						numFN--;
					}
				}

				precision = 1.0 * numTP / (numTP + numFP);
				recall = 1.0 * numTP / (numTP + numFN);
				if (!Double.isNaN(precision) && !Double.isNaN(recall)) {
					ret[i][0] = recall;
					ret[i][1] = precision;
				}
			}
		}
		log.info("PR curve created");
		return ret;

	}
}
