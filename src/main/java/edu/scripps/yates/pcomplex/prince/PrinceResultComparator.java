package edu.scripps.yates.pcomplex.prince;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.AbstractEpicOrPrinceResultsComparator;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.epic.EpicSummary;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.ProteinProteinInteraction;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.set.hash.THashSet;

public class PrinceResultComparator extends AbstractEpicOrPrinceResultsComparator {
	private final static Logger log = Logger.getLogger(PrinceResultComparator.class);
	public static final String clusterOne_results_file_suffix = "_interactions_all_clusterOne";
	private static List<ProteinComplexDB> dBs;

	public static void main(String[] args) {
		ProteinComplexAnalyzer.useComplexPortalDB = true;
		ProteinComplexAnalyzer.useCoreCorumDB = true;
		ProteinComplexAnalyzer.useHUMAP = true;

		final File princeClusterOneFile = new File(args[0]);
		final double minOverlapScore = Double.valueOf(args[1]);
		dBs = ProteinComplexAnalyzer.getDBs();
		final PrinceResultComparator epicComparator = new PrinceResultComparator(princeClusterOneFile, minOverlapScore);
		try {
			epicComparator.compareWithDBs(dBs);
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/**
	 * 
	 * @param princeClusterOneFile: file with the clusterOne results from prince.
	 *                              Use PrinceUtilities to read it.
	 * 
	 * @param minOverlapScore
	 */
	public PrinceResultComparator(File princeClusterOneFile, double minOverlapScore) {
		super(minOverlapScore, princeClusterOneFile);

	}

	/**
	 * Compare epic result against other EPIC result
	 * 
	 * @return
	 * 
	 * @throws IOException
	 */
	public File compareWithOtherPrinceResult(File princeClusterOneResultsFile) throws IOException {
		final File princeClusterOneFile = super.fileOrFolder;
		final String experimentName1 = getExperimentName();
		final List<ProteinComplex> complexes1 = readComplexes(princeClusterOneFile);
		Set<Object> set1 = complexes1.stream().map(c -> c.getId()).collect(Collectors.toSet());

		final String experimentName2 = getExperimentNameFromPrinceClusterOneResultsFile(princeClusterOneResultsFile);
		final List<ProteinComplex> complexes2 = readComplexes(princeClusterOneResultsFile);
		Set<Object> set2 = complexes2.stream().map(c -> c.getId()).collect(Collectors.toSet());

		final File outputFile = new File(
				princeClusterOneFile.getAbsolutePath() + File.separator + "epic_summary_" + minOverlapScore + ".txt");
		final FileWriter fw = new FileWriter(outputFile);

		// venn diagram for known Complexes
		set1 = complexes1.stream().filter(c -> c.isKnown()).map(c -> c.getId()).collect(Collectors.toSet());
		set2 = complexes2.stream().filter(c -> c.isKnown()).map(c -> c.getId()).collect(Collectors.toSet());
		final VennData vennKnown = new VennData("Known complexes", experimentName1, set1, experimentName2, set2, null,
				null);

		fw.write("Comparison between the KNOWN complexes between " + experimentName1 + " and " + experimentName2
				+ "\n\n");
		fw.write(vennKnown.toString() + "\n\n\n");

		// venn diagram for unknown Complexes
		final Set<ProteinComplex> unkown1 = complexes1.stream().filter(c -> !c.isKnown()).collect(Collectors.toSet());
		final Set<ProteinComplex> unkown2 = complexes2.stream().filter(c -> !c.isKnown()).collect(Collectors.toSet());
		final Set<ProteinComplex> overlapComplexes = new THashSet<ProteinComplex>();
		final Set<ProteinComplex> only1Complexes = new THashSet<ProteinComplex>();
		for (final ProteinComplex complex1 : unkown1) {
			boolean overlapped = false;
			for (final ProteinComplex complex2 : unkown2) {
				final double overlap = ClusterEvaluation.getOverlap(complex1, complex2);
				if (overlap >= minOverlapScore) {
					overlapped = true;
					break;
				}
			}
			if (overlapped) {
				// is in overlap
				overlapComplexes.add(complex1);
			} else {
				only1Complexes.add(complex1);
			}
		}
		final int only2Complexes = unkown2.size() - overlapComplexes.size();

		fw.write("Comparison between the UNKNOWN complexes between " + experimentName1 + " and " + experimentName2
				+ " with minimum overlap as " + minOverlapScore + "\n\n");
		final int union = overlapComplexes.size() + only1Complexes.size() + only2Complexes;
		fw.write("Union: " + union + " (" + getPercentageOverUnion(union, union) + "%)\n");
		fw.write("Unkown complexes from " + experimentName1 + ": " + unkown1.size() + " ("
				+ getPercentageOverUnion(unkown1.size(), union) + "%)\n");
		fw.write("Unkown complexes from " + experimentName2 + ": " + unkown2.size() + " ("
				+ getPercentageOverUnion(unkown2.size(), union) + "%)\n");
		fw.write("Unique to " + experimentName1 + ": " + only1Complexes.size() + " ("
				+ getPercentageOverUnion(only1Complexes.size(), union) + "%)\n");
		fw.write("Unique to " + experimentName2 + ": " + only2Complexes + " ("
				+ getPercentageOverUnion(only2Complexes, union) + "%)\n");
		fw.write("Overlap: " + overlapComplexes.size() + " (" + getPercentageOverUnion(overlapComplexes.size(), union)
				+ "%)\n");

		fw.write("Min overlap" + "\t" + "Overlapped complexes" + "\t" + "% Overlapped complexes" + "\n");
		double minOverlap = 0.0;
		while (minOverlap < 1.0) {
			try {

				final Set<ProteinComplex> overlapComplexesTMP = new THashSet<ProteinComplex>();

				for (final ProteinComplex complex1 : unkown1) {
					boolean overlapped = false;
					for (final ProteinComplex complex2 : unkown2) {
						final double overlap = ClusterEvaluation.getOverlap(complex1, complex2);
						if (overlap >= minOverlap) {
							overlapped = true;
							break;
						}
					}
					if (overlapped) {
						// is in overlap
						overlapComplexesTMP.add(complex1);
					}
				}
				fw.write(minOverlap + "\t" + overlapComplexesTMP.size() + "\t"
						+ getPercentageOverUnion(overlapComplexesTMP.size(), union) + "\n");
			} finally {
				minOverlap += 0.05;
			}
		}

		fw.write("\n\n\n\nComparison between all the complexes between " + experimentName1 + " and " + experimentName2
				+ " with minOverlap for unkown=" + minOverlapScore + "\n\n");
		final int unionTotal = vennKnown.getUnion123().size() + union;
		fw.write("Union: " + unionTotal + " (" + getPercentageOverUnion(unionTotal, unionTotal) + "%)\n");
		final int in1 = unkown1.size() + vennKnown.getSize1();
		fw.write("Unkown complexes from " + experimentName1 + ": " + in1 + " ("
				+ getPercentageOverUnion(in1, unionTotal) + "%)\n");
		final int in2 = unkown2.size() + vennKnown.getSize2();
		fw.write("Unkown complexes from " + experimentName2 + ": " + in2 + " ("
				+ getPercentageOverUnion(in2, unionTotal) + "%)\n");
		final int onlyIn1 = only1Complexes.size() + vennKnown.getUniqueTo1().size();
		fw.write("Unique to " + experimentName1 + ": " + onlyIn1 + " (" + getPercentageOverUnion(onlyIn1, unionTotal)
				+ "%)\n");
		final int onlyIn2 = only2Complexes + vennKnown.getUniqueTo2().size();
		fw.write("Unique to " + experimentName2 + ": " + onlyIn2 + " (" + getPercentageOverUnion(onlyIn2, unionTotal)
				+ "%)\n");
		final int overlap12 = overlapComplexes.size() + vennKnown.getIntersection12().size();
		fw.write("Overlap: " + overlap12 + " (" + getPercentageOverUnion(overlap12, unionTotal) + "%)\n");

		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
		return outputFile;
	}

	private static final DecimalFormat df = new DecimalFormat("#.#");

	private String getPercentageOverUnion(int num, int union) {
		return df.format(num * 100.0 / union);
	}

	/**
	 * Compare prince result against complex databases {@link ProteinComplexDB}
	 * 
	 * @throws IOException
	 */
	@Override
	public File compareWithDBs(List<ProteinComplexDB> dBs) throws IOException {
		final File princeClusterOneFile = fileOrFolder;
		final File outputFile = new File(
				princeClusterOneFile.getAbsolutePath() + File.separator + "prince_summary_" + minOverlapScore + ".txt");
		final FileWriter fw = new FileWriter(outputFile);

		writeSummaryHeaderLine(fw, dBs);

		final List<ProteinComplex> complexes = readComplexes(princeClusterOneFile);

		final String comparisonDescription = getComparisonDescription(complexes, null, dBs);
		fw.write(comparisonDescription);
		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
		return outputFile;

	}

	private void writeSummaryHeaderLine(FileWriter fw, List<ProteinComplexDB> dbs) throws IOException {
		fw.write("Experiment\t");
		fw.write("# complexes predicted\t");
		for (final ProteinComplexDB db : dbs) {
			fw.write(db.getName() + " # complexes known as " + minOverlapScore + "\t");
			fw.write(db.getName() + " # complexes unknown as " + minOverlapScore + "\t");
		}
		fw.write(EpicSummary.getHeaderLine());
		fw.write("\n");
	}

	@Override
	public String getExperimentName() {
		return getExperimentNameFromPrinceClusterOneResultsFile(fileOrFolder);
	}

	public String getExperimentNameFromPrinceClusterOneResultsFile(File clusterOneFileResults) {
		String name = FilenameUtils.getBaseName(clusterOneFileResults.getAbsolutePath());
		if (name.endsWith(clusterOne_results_file_suffix)) {
			name = name.substring(0, name.indexOf(clusterOne_results_file_suffix));
		}
		return name;
	}

	@Override
	public List<ProteinComplex> readComplexes(File complexesFile) throws IOException {
		final List<ProteinComplex> clusterOneResults = PrinceUtilities.readClusterOneResults(complexesFile);
		// now, evaluate if they are new or not by comparing with known complexes
		return clusterOneResults;
	}

	/**
	 * Reads protein-protein interactions from file typically called
	 * *._interactions_all.csv that is suppose to contain the protein-protein
	 * interactions from PRINCE results
	 * 
	 * @param expPredFile
	 * @return
	 * @throws IOException
	 */
	public List<ProteinProteinInteraction> readProteinProteinInteractions(File interactionsFile) throws IOException {
		final List<ProteinProteinInteraction> ret = new ArrayList<ProteinProteinInteraction>();

		final Stream<String> lines = Files.lines(interactionsFile.toPath());
		lines.forEachOrdered(line -> {
			final String[] split = line.split("\t");
			final String proteinA = split[2];
			final String proteinB = split[3];
			try {
				final ProteinComponent component1 = new ProteinComponent(proteinA, null);
				final ProteinComponent component2 = new ProteinComponent(proteinB, null);
				final double score = Double.valueOf(split[4]);
				final ProteinProteinInteraction ppi = new ProteinProteinInteraction(component1, component2, score);
				ret.add(ppi);
			} catch (final IOException e) {

			}
		});
		lines.close();
		log.info(ret.size() + " protein-protein interactions read from file " + interactionsFile.getAbsolutePath());
		return ret;
	}

}
