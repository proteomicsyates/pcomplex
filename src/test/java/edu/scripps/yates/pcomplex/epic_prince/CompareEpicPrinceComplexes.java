package edu.scripps.yates.pcomplex.epic_prince;

import java.io.File;
import java.io.IOException;
import java.util.List;

import edu.scripps.yates.pcomplex.AbstractEpicOrPrinceResultsComparator;
import edu.scripps.yates.pcomplex.epic.EpicResultComparator;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.prince.PrinceUtilities;

public class CompareEpicPrinceComplexes {
	private final File epicComplexFile;
	private final File princeComplexFile;
	private final File referenceComplexesFile;
	private final double minOverlapScore;

	/**
	 * First argument: file with epic complexes, which is a file named as
	 * Out.rf.*.clust.txt<br>
	 * Second argument: file with prince complexes, which is the output of
	 * clusterOne as a file named as: *_clusterOne_prec0.1.txt<br>
	 * Third argument: file with the reference complexes, that are collected by
	 * EPIC, as a file as: Out.rf.ref_complexes.txt<br>
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		final File epicComplexFile = new File(args[0]);
		final File princeComplexFile = new File(args[1]);
		final File referenceComplexesFile = new File(args[2]);
		final double minOverlapScore = 0.1;
		final CompareEpicPrinceComplexes c = new CompareEpicPrinceComplexes(epicComplexFile, princeComplexFile,
				referenceComplexesFile, minOverlapScore);
		try {
			c.run();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	public CompareEpicPrinceComplexes(File epicComplexFile, File princeComplexFile, File referenceComplexesFile,
			double minOverlapScore) {
		this.epicComplexFile = epicComplexFile;
		this.princeComplexFile = princeComplexFile;
		this.referenceComplexesFile = referenceComplexesFile;
		this.minOverlapScore = minOverlapScore;
	}

	public void run() throws IOException {
		final List<ProteinComplex> referenceComplexes = EpicResultComparator
				.readComplexesFromRefComplexesFile(referenceComplexesFile);
		final List<ProteinComplex> epicComplexes = EpicResultComparator.readComplexesFromClustFile(epicComplexFile);
		final List<ProteinComplex> princeComplexes = PrinceUtilities.readClusterOneResults(princeComplexFile);

		String description = AbstractEpicOrPrinceResultsComparator.getComparisonDescription(epicComplexes,
				princeComplexes, minOverlapScore);
		System.out.println("EPIC vs Prince:\n" + description);
		description = AbstractEpicOrPrinceResultsComparator.getComparisonDescription(epicComplexes, referenceComplexes,
				minOverlapScore);
		System.out.println("EPIC vs REF:\n" + description);
		description = AbstractEpicOrPrinceResultsComparator.getComparisonDescription(princeComplexes,
				referenceComplexes, minOverlapScore);
		System.out.println("PRINCE vs REF:\n" + description);
	}
}
