package edu.scripps.yates.pcomplex.epic;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;

public class EpicResultComparator {
	protected static final String suffix = "_out";
	private final File rootFolder;
	private int complexID = 0;
	private final double minOverlapScore;
	private static List<ProteinComplexDB> dBs;

	public EpicResultComparator(File file, double minOverlapScore) {
		rootFolder = file;
		this.minOverlapScore = minOverlapScore;
	}

	public static void main(String[] args) {
		ProteinComplexAnalyzer.useComplexPortalDB = true;
		ProteinComplexAnalyzer.useCoreCorumDB = true;
		ProteinComplexAnalyzer.useHUMAP = true;
		dBs = ProteinComplexAnalyzer.getDBs();
		final File folder = new File(args[0]);
		final double minOverlapScore = Double.valueOf(args[1]);
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			epicComparator.run();
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void run() throws IOException {

		if (rootFolder.isFile()) {
			throw new IllegalArgumentException(rootFolder.getAbsolutePath() + " is a file, not a folder");
		}
		final FileWriter fw = new FileWriter(
				new File(rootFolder.getAbsolutePath() + File.separator + "epic_summary_" + minOverlapScore + ".txt"));
		final File[] folders = getfilesEndingWith(rootFolder, suffix);
		writeSummaryHeaderLine(fw, dBs);
		for (final File folder : folders) {
			final String experimentName = getExperimentName(folder);
			final File complexesFile = getComplexesFile(folder);
			final List<ProteinComplex> complexes = readComplexes(complexesFile);

			final EpicSummary epicSummary = new EpicSummary(getSummaryFile(folder));
			writeSummaryLine(fw, experimentName, complexes, epicSummary, dBs);
		}
		fw.close();
	}

	private void writeSummaryLine(FileWriter fw, String experimentName, List<ProteinComplex> complexes,
			EpicSummary epicSummary, List<ProteinComplexDB> dBs) throws IOException {
		fw.write(experimentName + "\t");
		fw.write(complexes.size() + "\t");
		for (final ProteinComplexDB db : dBs) {
			int overlapping = 0;
			for (final ProteinComplex complex : complexes) {
				final double maxOverlap = complex.getMaxOverlap(db);
				if (maxOverlap > minOverlapScore) {
					overlapping++;
				}
			}
			fw.write(overlapping + "\t");
			fw.write((complexes.size() - overlapping) + "\t");
		}
		fw.write(epicSummary.getString());
		fw.write("\n");

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

	private String getExperimentName(File folder) {
		String name = FilenameUtils.getBaseName(folder.getAbsolutePath());
		if (name.endsWith(suffix)) {
			name = name.substring(0, name.indexOf(suffix));
		}
		return name;
	}

	private List<ProteinComplex> readComplexes(File complexesFile) throws IOException {
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

	private File[] getfilesEndingWith(File folder, String suffix) {
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

	private File getSummaryFile(File folder) {
		return getfilesEndingWith(folder, "Summary.txt")[0];
	}

	private File getComplexesFile(File folder) {
		return getfilesEndingWith(folder, "clust.txt")[0];
	}

}
