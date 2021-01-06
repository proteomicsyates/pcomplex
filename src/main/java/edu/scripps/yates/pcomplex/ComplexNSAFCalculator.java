package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.epic.EpicResultComparator;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Having a {@link SeparationExperiment} and a list of complexes (they may have
 * been derived by the analysis of the experiment by other software such as EPIC
 * or PRINCE), it will calculate the NSAF values of all the complexes
 * 
 * @author salvador
 *
 */
public class ComplexNSAFCalculator {

	private static final String uniprotFolder = "C:\\Users\\salvador\\Desktop\\uniprotKB";
	private static final String experimentFilesFolder = "C:\\Users\\salvador\\Desktop\\Anthony\\protein_complexes\\experiments";
	private final File complexListFile;
	private final SeparationExperiment separationExperiment;
	private File outputFile;
	private final UniprotProteinLocalRetriever uplr;
	private final ComplexFileType complexFileType;
	private boolean uniprotAnnotationsLoaded = false;
	private final TIntDoubleMap denominatorSumByFraction = new TIntDoubleHashMap();

	public static void main(String[] args) {
		try {
			// this can be the epic folder
			final File complexListFile = new File(args[0]);
			final SeparationExperiment separationExperiment = getSeparationExperiment(args[1]);
			final ComplexFileType complexFileType = ComplexFileType.valueOf(args[2]);
			final ComplexNSAFCalculator calculator = new ComplexNSAFCalculator(complexListFile, separationExperiment,
					complexFileType, uniprotFolder);

			final QuantType abundanceType = QuantType.NSAF;
			boolean printProteinsInDifferentLines = true;
			calculator.run(abundanceType, printProteinsInDifferentLines);
			printProteinsInDifferentLines = false;

			calculator.run(abundanceType, printProteinsInDifferentLines);
			System.out.println("Everything OK");
			System.exit(0);
		} catch (final IOException e) {

			e.printStackTrace();
			System.exit(-1);
		}
	}

	private static SeparationExperiment getSeparationExperiment(String experimentName) throws IOException {
		final SeparationExperiment ret = ProteinComplexAnalyzer.loadProjectSummaryFileNEW(experimentName,
				new File(experimentFilesFolder + File.separator + experimentName + ".tsv"));
		return ret;

	}

	public ComplexNSAFCalculator(File complexListFile, SeparationExperiment separationExperiment,
			ComplexFileType complexFileType, String uniprotFolder) {
		this.complexListFile = complexListFile;
		this.complexFileType = complexFileType;
		this.separationExperiment = separationExperiment;
		this.uplr = new UniprotProteinLocalRetriever(new File(uniprotFolder), true);
	}

	public void run(QuantType abundanceType, boolean printProteinsInDifferentLines) throws IOException {
		String pathname = null;
		if (complexListFile.isFile()) {
			pathname = complexListFile.getParent() + File.separator
					+ FilenameUtils.getBaseName(complexListFile.getAbsolutePath()) + "_sortedBy" + abundanceType;
		} else {
			pathname = complexListFile.getAbsolutePath() + File.separator
					+ FilenameUtils.getBaseName(complexListFile.getAbsolutePath()) + "_sortedBy" + abundanceType;
		}
		if (printProteinsInDifferentLines) {
			pathname += "_separated";
		}
		pathname += ".txt";
		this.outputFile = new File(pathname);
		final List<ProteinComplex> complexes = readProteinComplexes(this.complexListFile, complexFileType);
		for (final ProteinComplex proteinComplex : complexes) {
			final double abundance = calculateComplexAbundance(proteinComplex, separationExperiment, abundanceType);
			proteinComplex.setAbundance(abundance);
		}
		// sort and print sorted by abundance
		Collections.sort(complexes, new Comparator<ProteinComplex>() {

			@Override
			public int compare(ProteinComplex o1, ProteinComplex o2) {
				return Double.compare(o2.getAbundance(), o1.getAbundance());
			}
		});
		// print to file

		printToFile(complexes, outputFile, printProteinsInDifferentLines);
		System.out.println("File printed at " + outputFile.getAbsolutePath());
	}

	private List<ProteinComplex> readProteinComplexes(File epicFolder, ComplexFileType complexFileType)
			throws IOException {
		switch (complexFileType) {
		case EPIC:
			return EpicResultComparator.readPredictedProteinComplexes(epicFolder);

		case PRINCE:
			throw new IllegalArgumentException(complexFileType + " not supported...yet.");

		default:
			throw new IllegalArgumentException(complexFileType + " not known.");

		}
	}

	private void printToFile(List<ProteinComplex> complexes, File file, boolean printProteinsInDifferentLines2)
			throws IOException {
		if (printProteinsInDifferentLines2) {
			printToFileInSeparateLine(complexes, file);
		} else {
			printToFileInSameLine(complexes, file);
		}

	}

	private void printToFileInSameLine(List<ProteinComplex> complexes, File file) throws IOException {
		final FileWriter fw = new FileWriter(file);
		int rank = 1;
		for (final ProteinComplex proteinComplex : complexes) {
			fw.write(rank++ + "\t" + proteinComplex.getAbundance() + "\t" + proteinComplex.getId());
			final List<String> componentListKey = proteinComplex.getComponentListKey();
			for (final String accOrGene : componentListKey) {
				fw.write("\t" + accOrGene);
			}
			fw.write("\n");
		}
		fw.close();

	}

	private void printToFileInSeparateLine(List<ProteinComplex> complexes, File file) throws IOException {
		final FileWriter fw = new FileWriter(file);
		int rank = 1;
		for (final ProteinComplex proteinComplex : complexes) {

			final List<String> componentListKey = proteinComplex.getComponentListKey();
			for (final String accOrGene : componentListKey) {

				fw.write(rank + "\t" + proteinComplex.getAbundance() + "\t" + accOrGene + "\n");

			}
			rank++;
		}
		fw.close();

	}

	private double calculateComplexAbundance(ProteinComplex proteinComplex, SeparationExperiment separationExperiment,
			QuantType abundanceType) {
		switch (abundanceType) {
		case NSAF:
			return calculateComplexNSAF(proteinComplex, separationExperiment);

		case SPC:

			return calculateComplexSPC(proteinComplex, separationExperiment);
		case MW:
			return calculateComplexMW(proteinComplex, separationExperiment);
		default:
			throw new IllegalArgumentException("QuantType '" + abundanceType + "' not supported");
		}
	}

	private double calculateComplexMW(ProteinComplex proteinComplex, SeparationExperiment separationExperiment) {
		double mw = 0.0;
		for (final ProteinComponent component : proteinComplex.getComponentList()) {
			final List<Protein> proteinList = separationExperiment.getTotalProteinsByAcc().get(component.getKey());
			if (proteinList != null) {
				for (final Protein protein : proteinList) {
					if (protein.getMw() != null) {
						mw += protein.getMw();
						break;
					}
				}
			}
		}
		return mw;
	}

	private double calculateComplexSPC(ProteinComplex proteinComplex, SeparationExperiment separationExperiment) {
		double spc = 0.0;
		for (final ProteinComponent component : proteinComplex.getComponentList()) {
			final List<Protein> proteinList = separationExperiment.getTotalProteinsByAcc().get(component.getKey());
			if (proteinList != null) {
				for (final Protein protein : proteinList) {
					spc += protein.getSpc();
				}
			}
		}
		return spc;
	}

	/**
	 * NSAF for a complex is the sum of the NSAF of its components divided by the
	 * number of components in the complex
	 * 
	 * @param proteinComplex
	 * @param separationExperiment
	 * @return
	 */
	public double calculateComplexNSAF(ProteinComplex proteinComplex, SeparationExperiment separationExperiment) {
		double complexNSAF = 0.0;
		for (final ProteinComponent component : proteinComplex.getComponentList()) {
			final List<Fraction> fractions = separationExperiment.getFractions();
			for (final Fraction fraction : fractions) {
				if (fraction.containsReplicates()) {
					final TDoubleArrayList proteinNSAFPerReplicates = new TDoubleArrayList();
					for (final int replicate : fraction.getReplicates().toArray()) {
						final double nsafInReplicate = calculateProteinNSAF(component, fraction, replicate);
						if (!Double.isNaN(nsafInReplicate)) {
							proteinNSAFPerReplicates.add(nsafInReplicate);
						}
					}
					if (!proteinNSAFPerReplicates.isEmpty()) {
						// average
						complexNSAF += Maths.mean(proteinNSAFPerReplicates);
					}
				} else {
					final double proteinNSAF = calculateProteinNSAF(component, fraction);
					if (!Double.isNaN(proteinNSAF)) {
						complexNSAF += proteinNSAF;
					}
				}
			}
		}
		complexNSAF = complexNSAF / proteinComplex.getComponentList().size();
		return complexNSAF;
	}

	private double calculateProteinNSAF(ProteinComponent component, Fraction fraction) {
		return calculateProteinNSAF(component, fraction, 1);
	}

	/**
	 * NSAF is calculated as the number of spectral counts (SpC) identifying a
	 * protein, divided by the protein's length (L), divided by the sum of SpC/L for
	 * all proteins in the experiment.
	 * 
	 * @param component
	 * @param separationExperiment2
	 * @return
	 */
	private double calculateProteinNSAF(ProteinComponent component, Fraction fraction, int replicate) {
		// get annotations at once first
		if (!uniprotAnnotationsLoaded) {
			uplr.getAnnotatedProteins(null, separationExperiment.getProteinACCList());
			uniprotAnnotationsLoaded = true;
		}
		final double denominatorSum = getDenominatorSum(fraction, replicate);
		final String acc = component.getAcc();
		Integer spc = 0;

		final Protein protein = fraction.getProteinByAcc(acc);

		if (protein != null) {
			spc = protein.getSpc(replicate);
			if (spc == null) {
				return Double.NaN;
			}
		} else {
			return Double.NaN;
		}
		final double length = getLength(acc);
		final double numerator = spc * 1.0 / length;

		return numerator / denominatorSum;

	}

	public double calculateProteinNSAF(Protein protein, Fraction fraction, int replicate) {
		// get annotations at once first
		if (!uniprotAnnotationsLoaded) {
			uplr.getAnnotatedProteins(null, separationExperiment.getProteinACCList());
			uniprotAnnotationsLoaded = true;
		}
		final double denominatorSum = getDenominatorSum(fraction, replicate);

		Integer spc = 0;

		if (protein != null) {
			spc = protein.getSpc(replicate);
			if (spc == null) {
				return Double.NaN;
			}
		} else {
			return Double.NaN;
		}
		final String acc = protein.getAcc();
		final double length = getLength(acc);
		final double numerator = spc * 1.0 / length;

		return numerator / denominatorSum;

	}

	private double getDenominatorSum(Fraction fraction, int replicate) {

		if (!denominatorSumByFraction.containsKey(fraction.getFractionNumber())) {
			double denominatorSum = 0.0;
			for (final Protein protein : fraction.getProteins()) {
				final String acc = protein.getAcc();

				int spc = 0;

				if (protein != null) {
					spc += protein.getSpc(replicate);
				}

				final double length = getLength(acc);
				final double num = spc * 1.0 / length;

				if (!Double.isNaN(num)) {
					denominatorSum += num;
				}
			}
			denominatorSumByFraction.put(fraction.getFractionNumber(), denominatorSum);
		}
		return denominatorSumByFraction.get(fraction.getFractionNumber());
	}

	private double getLength(String acc) {
		final Set<String> accs = new THashSet<String>();
		if (acc.contains("#")) {
			final String[] split = acc.split("#");
			for (final String acc2 : split) {
				accs.add(acc2);
			}
		} else {
			accs.add(acc);
		}
		final TIntList lengths = new TIntArrayList();
		for (final String acc2 : accs) {

			final Map<String, Entry> annotatedProtein = this.uplr.getAnnotatedProtein(null, acc2);
			if (annotatedProtein != null) {
				final String proteinSequence = UniprotEntryUtil.getProteinSequence(annotatedProtein.get(acc2));
				if (proteinSequence != null) {
					lengths.add(proteinSequence.length());
				}
			}
		}
		return Maths.mean(lengths.toArray());
	}

}
