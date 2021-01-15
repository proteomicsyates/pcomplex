package edu.scripps.yates.pcomplex.prince;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.pcomplex.ClusterOneInterface;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;

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

	public static File createClusters(File interactionsFile) throws IOException {

		final ClusterOneInterface clusterOne = new ClusterOneInterface();
		final List<ProteinPairInteraction> interactions = new ArrayList<ProteinPairInteraction>();
		final List<String> lines = Files.readAllLines(interactionsFile.toPath());
		for (int i = 1; i < lines.size(); i++) {
			final String line = lines.get(i);
			final String[] split = line.split(",");
			int index = 1;
			if (split.length == 7) {
				index = 2;
			}
			final String proteinA = split[index];
			final String proteinB = split[index + 1];

			final double[] probabilities = new double[1];
			try {
				probabilities[0] = Double.valueOf(split[index + 2]);
				final ProteinPairInteraction ppi = new ProteinPairInteraction(proteinA, proteinB, probabilities,
						ClassLabel.INTRA_COMPLEX);
				interactions.add(ppi);
			} catch (final Exception e) {
				System.out.println(line);
			}
		}

		final double clusterOneMinDensity = 0.4;
		final double clusterOnePenalty = 2.9;
		final List<ProteinComplex> result = clusterOne.runClusterOne(interactions, clusterOneMinDensity,
				clusterOnePenalty);

		final File clustersFile = new File(interactionsFile.getParent() + File.separator
				+ FilenameUtils.getBaseName(interactionsFile.getAbsolutePath()) + "_clusterOne.txt");
		final FileWriter fw = new FileWriter(clustersFile);
		for (final ProteinComplex proteinComplex : result) {
			fw.write(proteinComplex.toString() + "\n");
		}
		fw.close();
		return clustersFile;
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
				componentsRaw = componentsRaw.substring(0);
			}
			if (componentsRaw.startsWith("]")) {
				componentsRaw = componentsRaw.substring(0, componentsRaw.length() - 1);
			}
			// split by '-'
			final String[] split2 = componentsRaw.split("-");
			for (String element : split2) {
				// remove "
				if (element.startsWith("\"")) {
					element = element.substring(0);
				}
				if (element.startsWith("\"")) {
					element = element.substring(0, element.length() - 1);
				}
				final ProteinComponent component = new ProteinComponent(element, element);
				complex.addComponent(component);
			}
			ret.add(complex);
		}
		return ret;
	}
}
