package edu.scripps.yates.pcomplex.epic;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.AbstractEpicOrPrinceResultsComparator;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.ProteinProteinInteraction;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.set.hash.THashSet;

public class EpicResultComparator extends AbstractEpicOrPrinceResultsComparator {
	private final static Logger log = Logger.getLogger(EpicResultComparator.class);
	protected static final String _out_suffix = "_out";
	private static int complexID = 0;
	private static List<ProteinComplexDB> dBs;

	public static void main(String[] args) {
		ProteinComplexAnalyzer.useComplexPortalDB = true;
		ProteinComplexAnalyzer.useCoreCorumDB = true;
		ProteinComplexAnalyzer.useHUMAP = true;

		final File folder = new File(args[0]);
		final double minOverlapScore = Double.valueOf(args[1]);
		dBs = ProteinComplexAnalyzer.getDBs();
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
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
	 * @param epicResultsFolder
	 *                          <ul>
	 *                          <li>this can be the folder in which multiple results
	 *                          are stored (all the folders ending on _out) which in
	 *                          that case it will be suitable for running the
	 *                          comparison against a {@link ProteinComplexDB}</li>
	 *                          <li>but it could be just a folder with one EPIC
	 *                          result in which case it is suitable for the
	 *                          comparison against another individual EPIC
	 *                          result</li>
	 *                          </ul>
	 * 
	 * 
	 * @param minOverlapScore
	 */
	public EpicResultComparator(File epicResultsFolder, double minOverlapScore) {
		super(minOverlapScore, epicResultsFolder);
		if (epicResultsFolder.isFile()) {
			throw new IllegalArgumentException(epicResultsFolder.getAbsolutePath() + " is a file, not a folder");
		}
	}

	/**
	 * Compare epic result against other EPIC result
	 * 
	 * @return
	 * 
	 * @throws IOException
	 */
	public File compareWithOtherEPICResult(File epicResultsFolder2) throws IOException {
		final File epicResultsFolder = super.fileOrFolder;
		final String experimentName1 = getExperimentName();
		final List<ProteinComplex> complexes1 = readPredictedProteinComplexes(epicResultsFolder);
		Set<Object> set1 = complexes1.stream().map(c -> c.getId()).collect(Collectors.toSet());

		final String experimentName2 = getExperimentNameFromEpicFolder(epicResultsFolder2);
		final List<ProteinComplex> complexes2 = readPredictedProteinComplexes(epicResultsFolder2);
		Set<Object> set2 = complexes2.stream().map(c -> c.getId()).collect(Collectors.toSet());

		final File outputFile = new File(
				epicResultsFolder.getAbsolutePath() + File.separator + "epic_summary_" + minOverlapScore + ".txt");
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
	 * Compare epic result againts complex databases {@link ProteinComplexDB}
	 * 
	 * @throws IOException
	 */
	@Override
	public File compareWithDBs(List<ProteinComplexDB> dBs) throws IOException {
		final File epicResultsFolder = super.fileOrFolder;
		if (epicResultsFolder.isFile()) {
			throw new IllegalArgumentException(epicResultsFolder.getAbsolutePath() + " is a file, not a folder");
		}
		final File outputFile = new File(
				epicResultsFolder.getAbsolutePath() + File.separator + "epic_summary_" + minOverlapScore + ".txt");
		final FileWriter fw = new FileWriter(outputFile);
		final File[] folders = getfilesEndingWith(epicResultsFolder, _out_suffix);
		writeSummaryHeaderLine(fw, dBs);
		for (final File folder : folders) {
			final File complexesFile = getPredictedComplexesFile(folder);
			final List<ProteinComplex> complexes = readComplexes(complexesFile);

			final EpicSummary epicSummary = new EpicSummary(getSummaryFile(folder));
			final String comparisonDescription = getComparisonDescription(complexes, epicSummary.getString(), dBs);
			fw.write(comparisonDescription);
		}
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
		return getExperimentNameFromEpicFolder(fileOrFolder);
	}

	private String getExperimentNameFromEpicFolder(File epicFolder) {
		String name = FilenameUtils.getBaseName(epicFolder.getAbsolutePath());
		if (name.endsWith(_out_suffix)) {
			name = name.substring(0, name.indexOf(_out_suffix));
		}
		return name;
	}

	@Override
	public List<ProteinComplex> readComplexes(File complexesFile) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final List<String> lines = Files.readAllLines(complexesFile.toPath());
		for (final String line : lines) {
			final ProteinComplex complex = new ProteinComplex(String.valueOf(++complexID));
			final String[] split = line.split("\t");
			for (final String protein : split) {
				final ProteinComponent component = new ProteinComponent(protein, null);
				complex.addComponent(component);
			}
			ret.add(complex);
		}
		return ret;
	}

	/**
	 * Reads protein complexes from file typically called Out.rf.ref_complexes.txt
	 * that is suppose to contain known protein complexes that are used as
	 * reference, having the complex ID name at the first column and all the
	 * components in the second column, separated by ","
	 * 
	 * @param complexesFile
	 * @return
	 * @throws IOException
	 */
	public static List<ProteinComplex> readComplexesFromRefComplexesFile(File complexesFile) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final List<String> lines = Files.readAllLines(complexesFile.toPath());
		for (final String line : lines) {
			final String[] split = line.split("\t");
			final ProteinComplex complex = new ProteinComplex(split[0]);
			complex.setKnown(true);
			for (int i = 1; i < split.length; i++) {

				final String proteins = split[i];
				final String[] split2 = proteins.split(",");
				for (final String protein : split2) {
					final ProteinComponent component = new ProteinComponent(protein, null);
					complex.addComponent(component);
				}

			}
			ret.add(complex);
		}
		log.info(ret.size() + " known protein complexes read from file " + complexesFile.getAbsolutePath());
		return ret;
	}

	public static List<ProteinComplex> readPredictedProteinComplexes(File epicFolder) throws IOException {

		if (epicFolder.isFile()) {

			final String fileName = FilenameUtils.getName(epicFolder.getAbsolutePath());
			if (fileName.endsWith("ref_complexes.txt")) {
				final List<ProteinComplex> knownComplexes = EpicResultComparator
						.readComplexesFromRefComplexesFile(epicFolder);
				return knownComplexes;
			} else if (fileName.endsWith("clust.txt")) {
				final List<ProteinComplex> unknownComplexes = EpicResultComparator
						.readComplexesFromClustFile(epicFolder);
				return unknownComplexes;
			}
			throw new IllegalArgumentException(epicFolder.getAbsolutePath()
					+ " MUST be a folder in which EPIC results are stored (Out.rf.ref_complexes.txt and Out.rf.exp.clust.txt) or one of these two files");
		}

		final List<ProteinComplex> referenceComplexes = EpicResultComparator
				.readComplexesFromRefComplexesFile(getReferenceComplexFile(epicFolder));
		final List<ProteinComplex> predictedComplexes = EpicResultComparator
				.readComplexesFromClustFile(getPredictedComplexesFile(epicFolder));

		// we mark as known the ones that have an overlap of at least 0.25 with at least
		// one known complex from the reference
		int known = 0;
		for (final ProteinComplex predicted : predictedComplexes) {
			for (final ProteinComplex reference : referenceComplexes) {
				final double overlap = ClusterEvaluation.getOverlap(predicted, reference);
				if (overlap >= 0.25) {
					predicted.setKnown(true);
					predicted.setId(predicted.getId() + " (overlap " + overlap + " with " + reference.getId() + ")");
					known++;
					break;
				}
			}
		}
		log.info(predictedComplexes.size() + " complexes, (" + known + " known and "
				+ (predictedComplexes.size() - known) + " novel) from EPIC results in " + epicFolder.getAbsolutePath());
		return predictedComplexes;

	}

	public static File getReferenceComplexFile(File epicFolder) {
		final File[] files = getfilesEndingWith(epicFolder, "ref_complexes.txt");
		if (files.length > 0) {
			return files[0];
		}
		return null;
	}

	private static File[] getfilesEndingWith(File folder, String suffix) {
		final File[] files = folder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.endsWith(suffix)) {
					return true;
				}
				return false;
			}
		});
		return files;
	}

	private static File getSummaryFile(File folder) {
		final File[] files = getfilesEndingWith(folder, "Summary.txt");
		if (files.length > 0) {
			return files[0];
		}
		return files[0];
	}

	private static File getPredictedComplexesFile(File folder) {
		final File[] files = getfilesEndingWith(folder, "clust.txt");
		if (files.length > 0) {
			return files[0];
		}
		return files[0];
	}

	/**
	 * Reads protein complexes from file typically called Out.rf.exp.clust.txt or
	 * Out.rf.comb.clust.txt that is suppose to contain unknown protein complexes
	 * from EPIC results, each line being a complex in which its components are
	 * separated in different columns
	 * 
	 * @param clustFile
	 * @return
	 * @throws IOException
	 */
	public static List<ProteinComplex> readComplexesFromClustFile(File clustFile) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final List<String> lines = Files.readAllLines(clustFile.toPath());
		int num = 1;
		for (final String line : lines) {

			final ProteinComplex complex = new ProteinComplex("Predicted-" + num++);
			complex.setKnown(false);
			final String[] split = line.split("\t");
			for (int i = 0; i < split.length; i++) {

				final String protein = split[i];

				final ProteinComponent component = new ProteinComponent(protein, null);
				complex.addComponent(component);

			}
			ret.add(complex);
		}
		log.info(ret.size() + " predicted protein complexes read from file " + clustFile.getAbsolutePath());
		return ret;
	}

	/**
	 * Reads protein-protein interactions from file typically called
	 * Out.rf.exp.pred.txt that is suppose to contain the protein-protein
	 * interactions from EPIC results
	 * 
	 * @param expPredFile
	 * @return
	 * @throws IOException
	 */
	public static List<ProteinProteinInteraction> readProteinProteinInteractionsFromExpClustFile(File expPredFile)
			throws IOException {
		final List<ProteinProteinInteraction> ret = new ArrayList<ProteinProteinInteraction>();
		final List<String> lines = Files.readAllLines(expPredFile.toPath());

		for (final String line : lines) {
			final String[] split = line.split("\t");
			final ProteinComponent component1 = new ProteinComponent(split[0], null);
			final ProteinComponent component2 = new ProteinComponent(split[1], null);
			final double score = Double.valueOf(split[2]);
			final ProteinProteinInteraction ppi = new ProteinProteinInteraction(component1, component2, score);
			ret.add(ppi);
		}
		log.info(ret.size() + " protein-protein interactions read from file " + expPredFile.getAbsolutePath());
		return ret;
	}

}
