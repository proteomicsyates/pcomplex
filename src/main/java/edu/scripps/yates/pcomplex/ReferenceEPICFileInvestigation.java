package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import edu.scripps.yates.pcomplex.epic.EpicResultComparator;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;

/**
 * We are not sure whether the file with the 'ref' results contain real
 * complexes in the sample or it is just a reference set use for training the ML
 * algorithm
 * 
 * @author salvador
 *
 */
public class ReferenceEPICFileInvestigation {

	private static File epicFolder = new File("C:\\Users\\salvador\\Desktop\\epic\\input\\Beta_cell_PCP_SPC_out");
	private static final double maxOverlapScore = 0.25;

	public static void main(String[] args) {
		List<ProteinComplex> complexes;
		try {
			complexes = EpicResultComparator.readPredictedProteinComplexes(epicFolder);

			final List<ProteinComplex> knownList = complexes.stream().filter(c -> c.isKnown())
					.collect(Collectors.toList());
			final List<ProteinComplex> unknownList = complexes.stream().filter(c -> !c.isKnown())
					.collect(Collectors.toList());

			final List<ProteinComplex> unknownFoundInReference = new ArrayList<ProteinComplex>();
			final List<ProteinComplex> knownFoundInUnknown = new ArrayList<ProteinComplex>();
			for (final ProteinComplex unknown : unknownList) {
				for (final ProteinComplex known : knownList) {
					final double overlap = ClusterEvaluation.getOverlap(unknown, known);
					if (overlap > maxOverlapScore) {
						unknownFoundInReference.add(unknown);
						knownFoundInUnknown.add(known);
					}
				}
			}
			for (int i = 0; i < unknownFoundInReference.size(); i++) {
				final ProteinComplex unknown = unknownFoundInReference.get(i);
				final ProteinComplex ref = knownFoundInUnknown.get(i);
				final double overlap = ClusterEvaluation.getOverlap(unknown, ref);
				System.out.println((i + 1) + "\tOverlap: " + overlap + "\nPredicted:\t" + unknown
						+ "\nKnown in reference:\t" + ref + "\n\n");
			}
			System.out.println(knownList.size() + " complexes in reference");
			System.out.println(unknownList.size() + " complexes predicted");
			System.out.println(unknownFoundInReference.size() + " complexes found in reference");

		} catch (final IOException e) {
			e.printStackTrace();
		}
	}
}
