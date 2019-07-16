package edu.scripps.yates.pcomplex.model;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class SeparationExperiment {

	private final List<Fraction> fractions = new ArrayList<Fraction>();
	private final Map<String, Fraction> fractionsByName = new THashMap<String, Fraction>();
	private final Map<String, Set<Fraction>> fractionNumbersByProteinAccs = new THashMap<String, Set<Fraction>>();
	private final String projectName;
	private final static Logger log = Logger.getLogger(SeparationExperiment.class);
	private final List<Protein> proteinList = new ArrayList<Protein>();
	private final List<String> proteinACCList = new ArrayList<String>();

	public SeparationExperiment(String projectName) {
		this.projectName = projectName;
	}

	public void addProtein(String fractionName, int fractionNumber, Protein protein) {
		proteinList.add(protein);
		if (fractionsByName.containsKey(fractionName)) {
			fractionsByName.get(fractionName).addProtein(protein);
		} else {
			final Fraction fraction = new Fraction(fractionName, fractionNumber);
			fraction.addProtein(protein);
			fractionsByName.put(fractionName, fraction);
			fractions.add(fraction);
		}
		if (fractionNumbersByProteinAccs.containsKey(protein.getAcc())) {
			fractionNumbersByProteinAccs.get(protein.getAcc()).add(fractionsByName.get(fractionName));
		} else {
			final Set<Fraction> fractions = new THashSet<Fraction>();
			fractions.add(fractionsByName.get(fractionName));
			fractionNumbersByProteinAccs.put(protein.getAcc(), fractions);
		}
	}

	public List<Fraction> getFractions() {
		return fractions;
	}

	public Map<String, Fraction> getFractionsByName() {
		return fractionsByName;
	}

	public String getProjectName() {
		return projectName;
	}

	public List<Fraction> getSortedFractions() {
		fractions.sort(new Comparator<Fraction>() {

			@Override
			public int compare(Fraction o1, Fraction o2) {
				return Integer.compare(o1.getFractionNumber(), o2.getFractionNumber());
			}
		});
		return fractions;
	}

	public List<Protein> getTotalProteins() {
		return proteinList;
	}

	/**
	 * Map with key as ACC and value as the list of proteins in each fraction in
	 * order
	 * 
	 * @param minFractions        minimum number of fractions in which a protein
	 *                            should be detected
	 * @param minSPCInOneFraction minimum number of SPCs to consider that a protein
	 *                            is present in a fraction
	 * @return
	 */
	public Map<String, List<Protein>> getTotalProteinsByAcc(int minFractions, int minSPCInOneFraction) {

		final List<Fraction> sortedFractions = getSortedFractions();
		final Map<String, List<Protein>> ret = new THashMap<String, List<Protein>>();

		final List<Protein> totalProteins = getTotalProteins();
		for (final Protein protein : totalProteins) {
			final String acc = protein.getAcc();
			if (!ret.containsKey(acc)) {
				final ArrayList<Protein> list = new ArrayList<Protein>(fractions.size());
				for (int i = 0; i < fractions.size(); i++) {
					list.add(null);
				}
				ret.put(acc, list);
			}
			final int indexOf = sortedFractions.indexOf(fractionsByName.get(protein.getFractionName()));
			final List<Protein> list = ret.get(acc);
			list.set(indexOf, protein);
		}
		if (minFractions > 0) {
			final List<String> toRemove = new ArrayList<String>();
			// filter by minSPC
			for (final String acc : ret.keySet()) {
				final List<Protein> list = ret.get(acc);
				int numFractions = 0;
				for (final Protein p : list) {
					if (p != null) {
						if (p.getSpc() >= minSPCInOneFraction) {
							numFractions++;
						}
					}
				}
				if (numFractions < minFractions) {
					toRemove.add(acc);
				}
			}
			if (toRemove.size() > 0) {
				for (final String acc : toRemove) {
					ret.remove(acc);
				}
				log.info("Number of proteins removed that were not in a min num fractions = " + minFractions
						+ " with spc>= " + minSPCInOneFraction + " were " + toRemove.size());
			}
		}
		return ret;
	}

	public Set<ProteinComplex> getCompleteComplexes(ProteinComplexDB proteinComplexDB, int minNumComponentsInComplex) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		for (final Fraction fraction : fractions) {
			ret.addAll(fraction.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex));
		}
		return ret;
	}

	public Set<ProteinComplex> getCompleteComplexes(ProteinComplexDB proteinComplexDB,
			ProteinComplexExistenceCriteria existenceCriteria, int minNumComponentsInComplex) {
		final Set<ProteinComplex> ret = new THashSet<ProteinComplex>();
		for (final Fraction fraction : fractions) {
			ret.addAll(fraction.getCompleteComplexes(proteinComplexDB, existenceCriteria, minNumComponentsInComplex));
		}
		return ret;
	}

	/**
	 * Exports a text table separated by tabs in a file located in the same folder
	 * as the properties file and with the name of the project.
	 * 
	 * @return
	 * @throws IOException
	 */
	public File exportToCSV(File folder, DataType dataType) throws IOException {
		final List<Fraction> fractionsSorted = getFractions();
		// sort fractions by fraction number
		Collections.sort(fractionsSorted, new Comparator<Fraction>() {

			@Override
			public int compare(Fraction o1, Fraction o2) {
				return Integer.compare(o1.getFractionNumber(), o2.getFractionNumber());
			}
		});
		// take all acccessions in the experiment
		final List<String> accs = new ArrayList<String>();
		accs.addAll(fractionNumbersByProteinAccs.keySet());
		// sort them alphabetically
		Collections.sort(accs);
		// open file to write
		final File outputFile = new File(
				folder.getAbsolutePath() + File.separator + projectName + "_" + dataType + ".tsv");
		final FileWriter fw = new FileWriter(outputFile);
		// header
		final String separator = ",";
		fw.write("ACC" + separator + "Replicate");
		for (final Fraction fraction : fractionsSorted) {
			fw.write(separator + fraction.getFractionNumber());
		}
		fw.write("\n");
		for (final String acc : accs) {
			boolean proteinInfoWritten = false;
			final StringBuilder sb = new StringBuilder();
			for (final Fraction fraction : fractionsSorted) {
				final Protein protein = fraction.getProteinByAcc(acc);
				if (!proteinInfoWritten && protein != null) {
					sb.insert(0, protein.getAcc() + separator + "1" + separator);
					proteinInfoWritten = true;
				}
				if (protein != null) {
					sb.append(getData(protein, dataType));
				} else {
					sb.append("NaN");
				}
				sb.append(separator);
			}
			fw.write(sb.toString() + "\n");
		}
		fw.close();
		return outputFile;
	}

	private String getData(Protein protein, DataType dataType) {
		switch (dataType) {
		case NSAF:
			return String.valueOf(protein.getNSAF());

		case SPC:
			return String.valueOf(protein.getSpc());

		default:
			break;
		}
		throw new IllegalArgumentException(dataType.name() + " not recognized");
	}
}
