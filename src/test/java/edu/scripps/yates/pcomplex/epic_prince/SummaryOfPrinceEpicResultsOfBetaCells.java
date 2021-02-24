package edu.scripps.yates.pcomplex.epic_prince;

import java.io.File;
import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.prince.PrinceResultComparator;
import edu.scripps.yates.pcomplex.prince.PrinceUtilities;

public class SummaryOfPrinceEpicResultsOfBetaCells {
	private final static String princeResultsFolderPath = "C:\\Users\\salvador\\Dropbox (Scripps Research)\\beta_cells_PCP\\prince";
	// to get the gold standard, we need the epic results
	private final static String epicResultsFolderPath = "C:\\Users\\salvador\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic";
	private final static String experiment_name = "Beta_cell_PCP_NSAF";

	private final static File findClusterOneResultFile(File interactionsFile, String suffix) {
		final File ret = new File(getPathToClusterOneResultsFile(interactionsFile, suffix));
		if (ret.exists()) {
			return ret;
		}
		return null;
	}

	private final static String getPathToClusterOneResultsFile(File interactionsFile, String suffix) {
		final String string = princeResultsFolderPath + File.separator
				+ FilenameUtils.getBaseName(interactionsFile.getAbsolutePath())
				+ PrinceResultComparator.prince_clusterOne_results_file_suffix + suffix + ".txt";
		return string;
	}

	private final static File findEpicResultsFolderBySpecies(String species) {
		final File epicResultsFolder = new File(
				epicResultsFolderPath + File.separator + "Beta_cell_PCP_NSAF" + "_" + species + "_out");
		if (epicResultsFolder.exists()) {
			return epicResultsFolder;
		}
		return null;
	}

	private final static DecimalFormat df = new DecimalFormat("#.#");

	/**
	 * First argument is the species either Human, Mouse or Rat<br>
	 * Second argument [OPTIONAL] is the precision cutoff or NaN<br>
	 * Third argument [OPTIONAL] is the bayes probability cutoff or NaN<br>
	 * Forth argument [OPTIONAL] is the classifier NB (by default), RF, SVM, LR
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			final String species = args[0].trim();
			Double precisionCutOff = null;
			Double bayesScoreCutOff = null;
			if (args.length > 1) {
				precisionCutOff = Double.valueOf(args[1].trim());
			}
			if (args.length > 2) {
				bayesScoreCutOff = Double.valueOf(args[2].trim());
			}
			String classifier = "NB";
			if (args.length > 3) {
				classifier = args[3].trim();
			}

			final File profilesFile = getPrinceProfilesFile(species);
			final File interactionsFile = findPrinceResultsInteractionsFile(profilesFile, classifier);
			final List<ProteinPairInteraction> ppis = PrinceUtilities.getProteinProteinInteractions(interactionsFile,
					precisionCutOff);
			final Double correspondingScoreThreshold = ppis.get(ppis.size() - 1).getProbabilities()[0];

			final String suffix = "_prec" + df.format(precisionCutOff);
			final File princeClusterOneResultFile = findClusterOneResultFile(interactionsFile, suffix);
			List<ProteinComplex> proteinComplexes = null;
			if (princeClusterOneResultFile == null || !princeClusterOneResultFile.exists()) {
				proteinComplexes = PrinceUtilities.createClusters(interactionsFile, precisionCutOff, bayesScoreCutOff,
						new File(getPathToClusterOneResultsFile(interactionsFile, suffix)));
			} else {
				proteinComplexes = PrinceUtilities.readClusterOneResults(princeClusterOneResultFile);
			}

			System.out.println(species + ":");
			System.out.println("PRINCE:");
			System.out.println("Classifier: " + classifier);
			System.out.println("Precision cutoff: " + precisionCutOff);
			System.out.println("Score cutoff: " + bayesScoreCutOff);
			System.out.println(princeClusterOneResultFile);
			System.out.println("Num complexes: " + proteinComplexes.size());
			System.out.println("corresponding score threshold: " + correspondingScoreThreshold);
			final File epicResultsFolder = findEpicResultsFolderBySpecies(species);
			System.out.println(epicResultsFolder);

			System.out.println("Everything correct");
		} catch (

		final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private static File getPrinceProfilesFile(String species) {
		return new File(princeResultsFolderPath + File.separator + experiment_name + "_" + species + ".csv");
	}

	private static File findPrinceResultsInteractionsFile(File profilesInputFile, String suffix) {
		final String pathToPrinceResultsInteractionsFile = getPathToPrinceResultsInteractionsFile(profilesInputFile,
				suffix);
		final File ret = new File(pathToPrinceResultsInteractionsFile);
		if (ret.exists()) {
			return ret;
		}
		return null;
	}

	private final static String getPathToPrinceResultsInteractionsFile(File profilesFile, String suffix) {
		final String string = princeResultsFolderPath + File.separator
				+ FilenameUtils.getBaseName(profilesFile.getAbsolutePath())
				+ PrinceResultComparator.prince_interactions_results_file_suffix + "_" + suffix + ".csv";
		return string;
	}
}
