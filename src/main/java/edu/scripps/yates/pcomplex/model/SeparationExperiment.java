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
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import smile.netlib.NLMatrix;

public class SeparationExperiment {

	private final List<Fraction> fractions = new ArrayList<Fraction>();
	private final Map<String, Fraction> fractionsByName = new THashMap<String, Fraction>();
	private final Map<String, Set<Fraction>> fractionNumbersByProteinAccs = new THashMap<String, Set<Fraction>>();
	private final String projectName;
	private final static Logger log = Logger.getLogger(SeparationExperiment.class);
	private static final int MIN_CONSECUTIVE_NON_ZERO_FRACTIONS = 2;
	private final List<Protein> proteinList = new ArrayList<Protein>();
	private final List<String> proteinACCList = new ArrayList<String>();
	private NLMatrix doubleNormalizedProfiles;
	private NLMatrix normalizedProfiles;
	private final Map<String, Set<String>> fractionNamesPerProteinAcc = new THashMap<String, Set<String>>();
	private File experimentFile;
	private Map<String, List<Protein>> totalProteinsByAcc;

	public SeparationExperiment(String projectName) {
		this.projectName = projectName;
	}

	public File getFile() {
		return experimentFile;
	}

	public void addProtein(String fractionName, int fractionNumber, Protein protein) {
		final String acc = protein.getAcc();
		if (!proteinACCList.contains(acc)) {
			proteinACCList.add(acc);
		}
		proteinList.add(protein);
		if (fractionsByName.containsKey(fractionName)) {
			fractionsByName.get(fractionName).addProtein(protein);
		} else {
			final Fraction fraction = new Fraction(fractionName, fractionNumber);
			fraction.addProtein(protein);
			fractionsByName.put(fractionName, fraction);
			fractions.add(fraction);
		}
		if (fractionNumbersByProteinAccs.containsKey(acc)) {
			fractionNumbersByProteinAccs.get(acc).add(fractionsByName.get(fractionName));
		} else {
			final Set<Fraction> fractions = new THashSet<Fraction>();
			fractions.add(fractionsByName.get(fractionName));
			fractionNumbersByProteinAccs.put(acc, fractions);
		}
	}

	public int getNumProteins() {
		return proteinACCList.size();
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

	/**
	 * Returns the total proteins detected in all fractions. Note that multiple
	 * objects referring to the same protein can be found, since a {@link Protein}
	 * object is a protein in a certain {@link Fraction}
	 * 
	 * @return
	 */
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

	/**
	 * Map with key as ACC and value as the list of proteins in each fraction in
	 * order
	 * 
	 * @return
	 */
	public Map<String, List<Protein>> getTotalProteinsByAcc() {
		if (totalProteinsByAcc == null) {
			totalProteinsByAcc = new THashMap<String, List<Protein>>();
			totalProteinsByAcc.putAll(getTotalProteinsByAcc(0, 0));
		}
		return totalProteinsByAcc;
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
	 * @param b
	 * 
	 * @return
	 * @throws IOException
	 */
	public File exportToTextSeparatedValues(File folder, DataType dataType, String separator,
			boolean includeReplicateColumn, boolean normalizeValuesAsMinEquals1) throws IOException {
		final List<Fraction> fractionsSorted = getFractions();

		// sort fractions by fraction number
		Collections.sort(fractionsSorted, new Comparator<Fraction>() {

			@Override
			public int compare(Fraction o1, Fraction o2) {
				return Integer.compare(o1.getFractionNumber(), o2.getFractionNumber());
			}
		});

		double multiplicativeFactor = 1.0;
		if (normalizeValuesAsMinEquals1) {
			// get the minimum non zero value and calculate the factor by which
			// we have to multiply
			double min = Double.MAX_VALUE;
			for (final String acc : fractionNumbersByProteinAccs.keySet()) {
				double minOfTheProtein = Double.MAX_VALUE;
				int numConsecutiveNonZeroFractions = 0;
				for (final Fraction fraction : fractionsSorted) {
					final Protein protein = fraction.getProteinByAcc(acc);
					if (protein != null) {
						final double data = protein.getData(dataType);
						if (!Double.isNaN(data)) {
							numConsecutiveNonZeroFractions++;
							if (minOfTheProtein > data) {
								minOfTheProtein = data;
							}
						} else {
							numConsecutiveNonZeroFractions = 0;
						}
					}
				}
				if (numConsecutiveNonZeroFractions >= MIN_CONSECUTIVE_NON_ZERO_FRACTIONS) {
					if (min > minOfTheProtein) {
						min = minOfTheProtein;
					}
				}
			}
			if (min < 1.0) {
				multiplicativeFactor = 1.0 / min;
			}
		}

		// take all acccessions in the experiment
		final List<String> accs = new ArrayList<String>();
		accs.addAll(fractionNumbersByProteinAccs.keySet());
		// sort them alphabetically
		Collections.sort(accs);
		// open file to write
		String extension = ".txt";
		if (separator == "\t") {
			extension = ".tsv";
		} else if (separator == ",") {
			extension = ".csv";
		}
		final File outputFile = new File(
				folder.getAbsolutePath() + File.separator + projectName + "_" + dataType + extension);
		final FileWriter fw = new FileWriter(outputFile);
		// header

		fw.write("ACC");
		if (includeReplicateColumn) {
			fw.write(separator + "Replicate");
		}
		for (final Fraction fraction : fractionsSorted) {
			fw.write(separator + fraction.getFractionNumber());
		}
		fw.write("\n");
		int numProteinDiscarded = 0;
		for (final String acc : accs) {
			boolean proteinInfoWritten = false;
			final StringBuilder sb = new StringBuilder();
			int maxNumNonZeroConsecutiveFractions = 0;
			int nonZeroFraction = 0;
			for (final Fraction fraction : fractionsSorted) {
				final Protein protein = fraction.getProteinByAcc(acc);
				if (!proteinInfoWritten && protein != null) {
					if (includeReplicateColumn) {
						sb.insert(0, protein.getAcc() + separator + "1" + separator);
					} else {
						sb.insert(0, protein.getAcc() + separator);
					}
					proteinInfoWritten = true;
				}
				if (protein != null) {
					sb.append(getData(protein, dataType, multiplicativeFactor));
					nonZeroFraction++;
					if (nonZeroFraction > maxNumNonZeroConsecutiveFractions) {
						maxNumNonZeroConsecutiveFractions = nonZeroFraction;
					}
				} else {
					sb.append("NaN");
					nonZeroFraction = 0;
				}
				sb.append(separator);
			}
			if (maxNumNonZeroConsecutiveFractions >= MIN_CONSECUTIVE_NON_ZERO_FRACTIONS) {
				fw.write(sb.toString() + "\n");
			} else {
				numProteinDiscarded++;
			}
		}
		fw.close();
		log.info(numProteinDiscarded + " proteins discarded for not having at least 2 consecutive non zero fractions");
		return outputFile;
	}

	private String getData(Protein protein, DataType dataType, double multiplicativeFactor) {
		switch (dataType) {
		case NSAF:
			return String.valueOf(protein.getNSAF() * multiplicativeFactor);

		case SPC:
			return String.valueOf(protein.getSpc() * multiplicativeFactor);

		default:
			break;
		}
		throw new IllegalArgumentException(dataType.name() + " not recognized");
	}

	public List<String> getProteinACCList() {
		return proteinACCList;
	}

	/**
	 * It gets the normalized data profile for protein acc1 in the
	 * {@link SeparationExperiment} exp.<br>
	 * If the experiment has not been normalized yet, it will be normalized in 2
	 * steps:<br>
	 * First, by normalizing each protein across all fractions (horizontally),<br>
	 * Second, by normalizing each fraction across all proteins (vertically).
	 * 
	 * @param exp
	 * @param acc1
	 * @param dataType
	 * @return
	 */
	public TDoubleList getNormalizedElutionProfile(String acc1, DataType dataType) {
		if (normalizedProfiles == null) {
			normalizedProfiles = new NLMatrix(getNumProteins(), getFractions().size(), 0.0);
			doubleNormalizedProfiles = new NLMatrix(getNumProteins(), getFractions().size(), 0.0);
			// populate matrix with the DataType of each protein in each
			// fraction
			for (int i = 0; i < getNumProteins(); i++) {
				final String acc = getProteinACCList().get(i);
				for (int j = 0; j < getFractions().size(); j++) {
					final Fraction fraction = getFractions().get(j);
					final Protein protein = fraction.getProteinByAcc(acc);
					if (protein != null) {
						normalizedProfiles.set(i, j, protein.getData(dataType));
					}
				}
			}

			// normalize by columns
			for (int j = 0; j < getFractions().size(); j++) {
				final TDoubleArrayList nums = new TDoubleArrayList();
				for (int i = 0; i < getNumProteins(); i++) {
					nums.add(normalizedProfiles.get(i, j));
				}
				for (int i = 0; i < getNumProteins(); i++) {
					normalizedProfiles.set(i, j, normalizedProfiles.get(i, j) / nums.sum());
					doubleNormalizedProfiles.set(i, j, normalizedProfiles.get(i, j) / nums.sum());
				}
			}
			// now normalize by rows
			for (int i = 0; i < getNumProteins(); i++) {
				final TDoubleArrayList nums = new TDoubleArrayList();
				for (int j = 0; j < getFractions().size(); j++) {
					nums.add(doubleNormalizedProfiles.get(i, j));
				}
				for (int j = 0; j < getFractions().size(); j++) {
					doubleNormalizedProfiles.set(i, j, doubleNormalizedProfiles.get(i, j) / nums.sum());
				}
			}
		}
		final TDoubleList ret = new TDoubleArrayList();
		final int i = getIndexForProtein(acc1);
		for (int j = 0; j < getFractions().size(); j++) {
			ret.add(normalizedProfiles.get(i, j));
		}
		return ret;
	}

	public TDoubleList getDoubleNormalizedElutionProfile(String acc1, DataType dataType) {
		getNormalizedElutionProfile(acc1, dataType);
		final TDoubleList ret = new TDoubleArrayList();
		final int i = getIndexForProtein(acc1);
		for (int j = 0; j < getFractions().size(); j++) {
			ret.add(doubleNormalizedProfiles.get(i, j));
		}
		return ret;
	}

	public int getIndexForProtein(String acc1) {
		return getProteinACCList().indexOf(acc1);
	}

	public TDoubleArrayList getElutionProfile(String acc1, DataType dataType) {
		final TDoubleArrayList ret = new TDoubleArrayList();
		for (final Fraction fraction : getSortedFractions()) {
			final Protein protein = fraction.getProteinByAcc(acc1);
			if (protein != null) {
				ret.add(protein.getData(dataType));
			} else {
				ret.add(0.0);
			}
		}
		return ret;
	}

	public Set<String> getFractionsInWhichProteinIsPresent(String acc) {
		if (!fractionNamesPerProteinAcc.containsKey(acc)) {
			final Set<String> fractions = getTotalProteinsByAcc().get(acc).stream().filter(p -> p != null)
					.map(p -> p.getFractionName()).filter(f -> f != null).collect(Collectors.toSet());
			fractionNamesPerProteinAcc.put(acc, fractions);
		}
		return fractionNamesPerProteinAcc.get(acc);
	}

	public void setFile(File file) {
		experimentFile = file;
	}
}
