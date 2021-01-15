package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.IOException;
import java.util.List;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;

public abstract class AbstractEpicOrPrinceResultsComparator {
	protected final double minOverlapScore;
	protected final File fileOrFolder;

	public AbstractEpicOrPrinceResultsComparator(double minOverlapScore, File fileOrFolder) {
		this.minOverlapScore = minOverlapScore;
		this.fileOrFolder = fileOrFolder;
	}

	public abstract String getExperimentName();

	public abstract List<ProteinComplex> readComplexes(File complexesFile) throws IOException;

	/**
	 * Gets a line comparing n DBs of complexes with the following values :<br>
	 * experiment name - total number of complexes we compare against the DBs -
	 * [number of complexes overlapping - number of complexes not overlapping with
	 * DB1 to n] - summary last column
	 * 
	 * @param complexes
	 * @param summaryLastColumn
	 * @param dBs
	 * @return
	 * @throws IOException
	 */
	protected String getComparisonDescription(List<ProteinComplex> complexes, String summaryLastColumn,
			List<ProteinComplexDB> dBs) throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(getExperimentName() + "\t");
		sb.append(complexes.size() + "\t");
		for (final ProteinComplexDB db : dBs) {
			int overlapping = 0;
			for (final ProteinComplex complex : complexes) {
				final double maxOverlap = complex.getMaxOverlap(db);
				if (maxOverlap > minOverlapScore) {
					overlapping++;
				}
			}
			sb.append(overlapping + "\t");
			sb.append((complexes.size() - overlapping) + "\t");
		}
		if (summaryLastColumn != null) {
			sb.append(summaryLastColumn);
		}
		sb.append("\n");
		return sb.toString();
	}

	/**
	 * Gets a line comparing a list of other complexes with the following values
	 * :<br>
	 * total number of complexes we want to compare against others - number of
	 * complexes overlapping with other complexes - number of complexes not
	 * overlapping with other complexes
	 * 
	 * @param complexes
	 * @param otherComplexes
	 * @return
	 */
	protected String getComparisonDescription(List<ProteinComplex> complexes, List<ProteinComplex> otherComplexes) {
		final StringBuilder sb = new StringBuilder();

		sb.append(complexes.size() + "\t");

		int overlapping = 0;
		for (final ProteinComplex complex : complexes) {
			final double maxOverlap = complex.getMaxOverlap(otherComplexes);
			if (maxOverlap > minOverlapScore) {
				overlapping++;
			}
		}
		sb.append(overlapping + "\t");
		sb.append((complexes.size() - overlapping) + "\t");

		sb.append("\n");
		return sb.toString();
	}

	/**
	 * Compare epic or prince result againts complex databases
	 * {@link ProteinComplexDB}
	 * 
	 * @throws IOException
	 */
	public abstract File compareWithDBs(List<ProteinComplexDB> dBs) throws IOException;
}
