package edu.scripps.yates.pcomplex.prince;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import edu.scripps.yates.pcomplex.ClusterOneInterface;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class PrinceUtilities {
	private final File interactionsFile = new File(
			"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\prince\\Beta_cell_PCP_NSAF_Human_interactions_all.csv");

	/**
	 * Reads the interactions created by PrinCE and apply the ClusterONE algorithm
	 * with the default parameters 0.4 and 2.9 and writes an output file with suffix
	 * _clusterOne.txt
	 * 
	 * @throws IOException
	 */

	public static List<ProteinComplex> createClusters(File interactionsFile, Double precisionCutOff, Double scoreCutOff,
			File clusterOneResultsFile) throws IOException {
		if (precisionCutOff != null && scoreCutOff != null) {
			throw new IllegalArgumentException("Either precision or score cut offs should be provided");
		}

		final List<ProteinPairInteraction> interactions = new ArrayList<ProteinPairInteraction>();
		final List<String> lines = Files.readAllLines(interactionsFile.toPath());
		final TObjectIntMap<String> indexesByColumnNames = new TObjectIntHashMap<String>();
		final String headers = lines.get(0).replace("\"", "");
		final String[] splitHeader = headers.split(",");
		for (int index = 0; index < splitHeader.length; index++) {
			final String header = splitHeader[index];
			indexesByColumnNames.put(header, index);
		}
		// depending on whether Prince was run with ensemble classifier or not, we will
		// have mean or score as the column with the score. the score is the higher the
		// better.
		int scoreIndex;
		if (indexesByColumnNames.containsKey("score")) {
			scoreIndex = indexesByColumnNames.get("score");
		} else {
			scoreIndex = indexesByColumnNames.get("mean");
		}
		// first we figure out in which line we want to set the threshold depending on
		// whether is a threshold over the score or over the precision
		//
		int maxLine = 1;
		for (int i = 1; i < lines.size(); i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");
			if (precisionCutOff != null) {
				// it would be in the maximum file line in which the precision is above the
				// threshold
				// lines before could have a lower precision, or even no precision at all, but
				// the interaction score will be higher and so we will keep them
				if (split.length < indexesByColumnNames.get("precision")) {
					continue;
				}
				final String precisionString = split[indexesByColumnNames.get("precision")];
				final Double precision = Double.valueOf(precisionString);
				if (precision >= precisionCutOff) {
					maxLine = i;
				}
			} else if (scoreCutOff != null) {
				final String scoreString = split[scoreIndex];
				final double score = Double.valueOf(scoreString);
				if (score >= scoreCutOff) {
					maxLine = i;
				}
			} else {
				throw new IllegalArgumentException("Either precision or score should be provided");
			}
		}
		System.out.println(
				"we will keep " + (maxLine - 1) + " interactions at that precision threshold " + precisionCutOff);

		for (int i = 1; i <= maxLine; i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");

			final String proteinA = split[indexesByColumnNames.get("protein_A")];
			final String proteinB = split[indexesByColumnNames.get("protein_B")];

			final double[] scores = new double[1];
			try {

				scores[0] = Double.valueOf(split[scoreIndex]);
				final ProteinPairInteraction ppi = new ProteinPairInteraction(proteinA, proteinB, scores,
						ClassLabel.INTRA_COMPLEX);
				interactions.add(ppi);
			} catch (final Exception e) {
				System.out.println(line);
			}
		}

		final ClusterOneInterface clusterOne = new ClusterOneInterface();
		final double clusterOneMinDensity = 0.4;
		final double clusterOnePenalty = 2.9;
		final List<ProteinComplex> result = clusterOne.runClusterOne(interactions, clusterOneMinDensity,
				clusterOnePenalty);

		final FileWriter fw = new FileWriter(clusterOneResultsFile);
		for (final ProteinComplex proteinComplex : result) {
			fw.write(proteinComplex.toString() + "\n");
		}
		fw.close();
		return result;
	}

	/**
	 * Reads cluster one result file which has a line per complex with this
	 * format:<br>
	 * 273100247 : ["P01308"-"P30410"-"Q6YK33"-"Q8HXV2"]
	 * 
	 * @param clusterOneFile
	 * @return
	 * @throws IOException
	 */
	public static List<ProteinComplex> readClusterOneResults(File clusterOneFile) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final List<String> lines = Files.readAllLines(clusterOneFile.toPath());
		for (final String line : lines) {
			final String[] split = line.split(":");
			final String complexID = split[0].trim();
			final ProteinComplex complex = new ProteinComplex(complexID);
			String componentsRaw = split[1].trim();
			// remove brackets
			if (componentsRaw.startsWith("[")) {
				componentsRaw = componentsRaw.substring(1);
			}
			if (componentsRaw.endsWith("]")) {
				componentsRaw = componentsRaw.substring(0, componentsRaw.length() - 1);
			}
			// split by '-'
			final String[] split2 = componentsRaw.split("-");
			for (String element : split2) {
				// remove "
				if (element.startsWith("\"")) {
					element = element.substring(1);
				}
				if (element.endsWith("\"")) {
					element = element.substring(0, element.length() - 1);
				}
				final ProteinComponent component = new ProteinComponent(element, element);
				complex.addComponent(component);
			}
			ret.add(complex);
		}
		return ret;
	}

	/**
	 * Reads an interactions file from Prince and makes a precision cut off and
	 * returns the corresponding score at that position
	 * 
	 * 
	 * @param interactionsFile
	 * @param precisionCutOff
	 * @return
	 * @throws IOException
	 */
	public static double getCorrespondingScoreThreshold(File interactionsFile, double precisionCutOff)
			throws IOException {
		final List<String> lines = Files.readAllLines(interactionsFile.toPath());
		final TObjectIntMap<String> indexesByColumnNames = new TObjectIntHashMap<String>();
		final String headers = lines.get(0);
		final String[] splitHeader = headers.split(",");
		for (int index = 0; index < splitHeader.length; index++) {
			final String header = splitHeader[index].replace("\"", "");
			indexesByColumnNames.put(header, index);
		}
		// depending on whether Prince was run with ensemble classifier or not, we will
		// have mean or score as the column with the score. the score is the higher the
		// better.
		int scoreIndex;
		if (indexesByColumnNames.containsKey("score")) {
			scoreIndex = indexesByColumnNames.get("score");
		} else {
			scoreIndex = indexesByColumnNames.get("mean");
		}
		int maxLine = 1;
		for (int i = 1; i < lines.size(); i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");

			// it would be in the maximum file line in which the precision is above the
			// threshold
			// lines before could have a lower precision, or even no precision at all, but
			// the interaction score will be higher and so we will keep them
			if (split.length < indexesByColumnNames.get("precision")) {
				continue;
			}
			final String precisionString = split[indexesByColumnNames.get("precision")];
			final Double precision = Double.valueOf(precisionString);
			if (precision >= precisionCutOff) {
				maxLine = i;
			}
		}

		final String line = lines.get(maxLine);
		final String[] split = line.split(",");
		final String scoreString = split[scoreIndex];
		final double score = Double.valueOf(scoreString);
		return score;
	}

	public static List<ProteinPairInteraction> getProteinProteinInteractions(File interactionsFile,
			double precisionCutOff) throws IOException {
		final List<String> lines = Files.readAllLines(interactionsFile.toPath());
		final TObjectIntMap<String> indexesByColumnNames = new TObjectIntHashMap<String>();
		final String headers = lines.get(0);
		final String[] splitHeader = headers.split(",");
		for (int index = 0; index < splitHeader.length; index++) {
			final String header = splitHeader[index].replace("\"", "");
			indexesByColumnNames.put(header, index);
		}
		// depending on whether Prince was run with ensemble classifier or not, we will
		// have mean or score as the column with the score. the score is the higher the
		// better.
		int scoreIndex;
		if (indexesByColumnNames.containsKey("score")) {
			scoreIndex = indexesByColumnNames.get("score");
		} else {
			scoreIndex = indexesByColumnNames.get("mean");
		}
		int maxLine = 1;
		for (int i = 1; i < lines.size(); i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");

			// it would be in the maximum file line in which the precision is above the
			// threshold
			// lines before could have a lower precision, or even no precision at all, but
			// the interaction score will be higher and so we will keep them
			if (split.length < indexesByColumnNames.get("precision")) {
				continue;
			}
			final String precisionString = split[indexesByColumnNames.get("precision")];
			final Double precision = Double.valueOf(precisionString);
			if (precision >= precisionCutOff) {
				maxLine = i;
			}
		}
		final List<ProteinPairInteraction> interactions = new ArrayList<ProteinPairInteraction>();
		for (int i = 1; i <= maxLine; i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");
			final String proteinA = split[indexesByColumnNames.get("protein_A")];
			final String proteinB = split[indexesByColumnNames.get("protein_B")];

			final double[] scores = new double[1];
			try {

				scores[0] = Double.valueOf(split[scoreIndex]);
				final ProteinPairInteraction ppi = new ProteinPairInteraction(proteinA, proteinB, scores,
						ClassLabel.INTRA_COMPLEX);
				interactions.add(ppi);
			} catch (final Exception e) {
				System.out.println(line);
			}
		}
		return interactions;
	}
}
