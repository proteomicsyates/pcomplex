package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.epic.EpicResultComparator;
import edu.scripps.yates.pcomplex.gprofiler.GProfiler;
import edu.scripps.yates.pcomplex.gprofiler.GProfilerResult;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.THashSet;

public class FiguresGenerator {
	private final static Logger log = Logger.getLogger(FiguresGenerator.class);
	private static final String GO_CYTOPLASM = "GO:0005737";
	private static final String GO_ORGANELLE = "GO:0043226";

	public static void main(String[] args) {
		final File epicFolder = new File(args[0]);
		final String experimentFile = args[1];

		final FiguresGenerator g = new FiguresGenerator(epicFolder, experimentFile);
		try {
			g.run();
			log.info("Everything OK");
			System.exit(0);
		} catch (final IOException e) {

			e.printStackTrace();
			System.exit(-1);
		}

	}

	private SeparationExperiment separationExperiment;
	private final File epicFolder;
	private final String experimentName;
	private final File projectSummaryFile;

	private void run() throws IOException {
		final File excelFile = getExcelFile();
		if (excelFile.exists()) {
			final boolean deleted = excelFile.delete();
			if (!deleted) {
				throw new IllegalArgumentException("unable to delete excel file at " + excelFile.getAbsolutePath());
			}
			log.info(excelFile.getAbsolutePath() + " deleted");
		}
		generateFig1A();

////				generateFig1B();
		generateFig1C();
		generateFig1D();
		generateFigSize();
	}

	private void generateFigSize() throws IOException {

		final TIntList sizes = new TIntArrayList();
		final List<ProteinComplex> complexes = EpicResultComparator.readProteinComplexes(epicFolder);
		for (final ProteinComplex proteinComplex : complexes) {
			sizes.add(proteinComplex.getComponentList().size());
		}
		final double avgSize = Maths.mean(sizes);
		final int smallCutOff = Double.valueOf(Math.floor(avgSize)).intValue();
		final List<ProteinComplex> small = complexes.stream().filter(c -> c.getComponentList().size() <= smallCutOff)
				.collect(Collectors.toList());
		final List<ProteinComplex> big = complexes.stream().filter(c -> c.getComponentList().size() > smallCutOff)
				.collect(Collectors.toList());
		if (small.size() + big.size() != complexes.size()) {
			throw new IllegalArgumentException("Split by sizes > or <= than " + smallCutOff + " is wrong");
		}
		final Set<String> proteinsInSmall = new THashSet<String>();
		for (final ProteinComplex complex : small) {
			complex.getComponentList().stream().forEach(c -> proteinsInSmall.add(c.getAcc()));
		}
		final Set<String> proteinsInBig = new THashSet<String>();
		for (final ProteinComplex complex : big) {
			complex.getComponentList().stream().forEach(c -> proteinsInBig.add(c.getAcc()));
		}
		final VennData venn = new VennData("Num proteins in Big or Small complexes", "small", proteinsInSmall, "big",
				proteinsInBig, null, null);

		final File figFile = File.createTempFile("figSize", ".tmp");
		figFile.deleteOnExit();
		final FileWriter fw = new FileWriter(figFile);
		fw.write(complexes.size() + " total complexes\n");
		fw.write("Avg size of complexes is " + avgSize + "\n");
		fw.write(small.size() + " small complexes (<=" + smallCutOff + ")\n");
		fw.write(big.size() + " big complexes (>" + smallCutOff + ")\n\n\n");
		fw.write(venn.toString());

		fw.close();
		FileUtils.separatedValuesToXLSX(figFile.getAbsolutePath(), getExcelFile().getAbsolutePath(), "\t", "FigSize");

		// now write files of big and small proteins
		final String expName = FilenameUtils.getBaseName(getExcelFile().getParentFile().getAbsolutePath());
		final File smallFile = new File(getExcelFile().getParentFile() + File.separator + expName
				+ "_Proteins_Lt_or_Eq_" + smallCutOff + ".txt");
		final FileWriter fw2 = new FileWriter(smallFile);
		for (final String acc : proteinsInSmall) {
			if (acc.contains("#")) {
				final String[] split = acc.split("#");
				for (final String acc2 : split) {
					fw2.write(acc2 + "\n");
				}
			} else {
				fw2.write(acc + "\n");
			}
		}
		fw2.close();
		log.info("File with proteins in small complexes is at: " + smallFile.getAbsolutePath());

		final File bigFile = new File(
				getExcelFile().getParentFile() + File.separator + expName + "_Proteins_Gt_" + smallCutOff + ".txt");
		final FileWriter fw3 = new FileWriter(bigFile);
		for (final String acc : proteinsInBig) {
			if (acc.contains("#")) {
				final String[] split = acc.split("#");
				for (final String acc2 : split) {
					fw3.write(acc2 + "\n");
				}
			} else {
				fw3.write(acc + "\n");
			}

		}
		fw3.close();
		log.info("File with proteins in big complexes is at: " + bigFile.getAbsolutePath());

	}

	private void generateFig1D() throws IOException {
		final List<ProteinComplex> complexes = EpicResultComparator.readProteinComplexes(epicFolder);

		final List<ProteinComplex> known = complexes.stream().filter(c -> c.isKnown()).collect(Collectors.toList());
		final List<ProteinComplex> unknown = complexes.stream().filter(c -> !c.isKnown()).collect(Collectors.toList());
		final StringBuilder sb = new StringBuilder();
		sb.append("Number of known complexes\t" + known.size() + "\n");
		sb.append("Number of new/unknown complexes\t" + unknown.size() + "\n");
		getSeparationExperiment();
		final List<Fraction> fractions = separationExperiment.getSortedFractions();

		// header
		for (final Fraction fraction : fractions) {
			sb.append("\t" + fraction.getFractionNumber());
		}
		sb.append("\n");

		// number of protein known protein complexes in fraction
		sb.append("# known protein complexes");
		for (final Fraction fraction : fractions) {
			int num = 0;
			for (final ProteinComplex complex : known) {
				if (complex.isFullyRepresented(fraction.getProteinsAndGenes())) {
					num++;
				}
			}
			sb.append("\t" + num);
		}
		sb.append("\n");
		// number of protein unknown protein complexes in fraction
		sb.append("# unknown protein complexes");
		for (final Fraction fraction : fractions) {
			int num = 0;
			for (final ProteinComplex complex : unknown) {
				if (complex.isFullyRepresented(fraction.getProteinsAndGenes())) {
					num++;
				}
			}
			sb.append("\t" + num);
		}
		sb.append("\n");
		// avg number of components
		sb.append("# avg protein complex size (number of components)");
		for (final Fraction fraction : fractions) {
			final TIntList sizes = new TIntArrayList();
			for (final ProteinComplex complex : complexes) {
				if (complex.isFullyRepresented(fraction.getProteinsAndGenes())) {
					sizes.add(complex.getComponentList().size());
				}
			}
			double mean = 0.0;
			if (!sizes.isEmpty()) {
				mean = Maths.mean(sizes);
			}
			sb.append("\t" + mean);
		}
		sb.append("\n");
		// avg MW of components
		sb.append("# avg MW (kDa) of complexes ");
		for (final Fraction fraction : fractions) {
			final TDoubleList mws = new TDoubleArrayList();
			for (final ProteinComplex complex : complexes) {

				if (complex.isFullyRepresented(fraction.getProteinsAndGenes())) {
					double complexMW = 0.0;
					for (final ProteinComponent component : complex.getComponentList()) {
						final Double mw = fraction.getProteinByAcc(component.getAcc()).getMw();
						if (mw != null) {
							complexMW += mw;
						} else {
							log.info("Protein " + component.getAcc() + " has no mw");
						}
					}
					mws.add(complexMW);
				}

			}
			double mean = 0.0;
			if (!mws.isEmpty()) {
				mean = Maths.mean(mws);
			}
			sb.append("\t" + mean / 1000.0);
		}
		sb.append("\n");

		final File fig1DFile = File.createTempFile("fig1D", ".tmp");
		fig1DFile.deleteOnExit();
		final FileWriter fw = new FileWriter(fig1DFile);
		fw.write(sb.toString() + "\n\n\n");

		fw.close();
		FileUtils.separatedValuesToXLSX(fig1DFile.getAbsolutePath(), getExcelFile().getAbsolutePath(), "\t", "Fig1D");
	}

	private void generateFig1C() throws IOException {
		final List<ProteinComplex> complexes = EpicResultComparator.readProteinComplexes(epicFolder);
		final List<ProteinComplex> known = complexes.stream().filter(c -> c.isKnown()).collect(Collectors.toList());
		final List<ProteinComplex> unknown = complexes.stream().filter(c -> !c.isKnown()).collect(Collectors.toList());
		final StringBuilder sb = new StringBuilder();
		sb.append("Number of known complexes\t" + known.size() + "\n");
		sb.append("Number of new/unknown complexes\t" + unknown.size() + "\n");

		// veen of GO vs CORUM vs IntAct
		final List<ProteinComplex> go = known.stream().filter(f -> f.isSourceQuickGO()).collect(Collectors.toList());
		final List<ProteinComplex> corum = known.stream().filter(f -> f.isSourceCorum()).collect(Collectors.toList());
		final List<ProteinComplex> intAct = known.stream().filter(f -> f.isSourceIntAct()).collect(Collectors.toList());
		final VennData venn = new VennData("QuickGO vs CORUM vs IntAct", "QuickGO", go, "CORUM", corum, "IntAct",
				intAct);

		final File fig1CFile = File.createTempFile("fig1C", ".tmp");
		fig1CFile.deleteOnExit();
		final FileWriter fw = new FileWriter(fig1CFile);
		fw.write(sb.toString() + "\n\n\n");
		fw.write(venn.toString());
		fw.close();
		FileUtils.separatedValuesToXLSX(fig1CFile.getAbsolutePath(), getExcelFile().getAbsolutePath(), "\t", "Fig1C");
	}

	private void generateFig1A() throws IOException {
		// from epic folder, take
		final File epicTotalProteinsGOCCFile = getEpicTotalProteinsGOCCFile();
		final GProfilerResult gProfilerResults = GProfiler.readGProfilerResults(epicTotalProteinsGOCCFile);
		final Set<String> cytoplasm = gProfilerResults.getProteinsWithGOTerm(GO_CYTOPLASM);
		final Set<String> organelle = gProfilerResults.getProteinsWithGOTerm(GO_ORGANELLE);
		final VennData venn = new VennData("Cytoplasm vs Organelle", "Cytoplasm", cytoplasm, "Organelle", organelle,
				null, null);
		final File fig1File = File.createTempFile("fig1A", ".tmp");
		fig1File.deleteOnExit();
		final FileWriter fw = new FileWriter(fig1File);
		fw.write(venn.toString());
		fw.close();
		FileUtils.separatedValuesToXLSX(fig1File.getAbsolutePath(), getExcelFile().getAbsolutePath(), "\t", "Fig1A");

	}

	private void generateFig1B() throws IOException {
		// from epic folder, take

		final List<ProteinComplex> complexes = EpicResultComparator.readProteinComplexes(epicFolder);
		final List<ProteinComplex> known = complexes.stream().filter(c -> c.isKnown()).collect(Collectors.toList());
		final List<ProteinComplex> unknown = complexes.stream().filter(c -> !c.isKnown()).collect(Collectors.toList());

		final File epicTotalProteinsGOCCFile = getEpicTotalProteinsGOCCFile();
		final GProfilerResult gProfilerResults = GProfiler.readGProfilerResults(epicTotalProteinsGOCCFile);

		final Set<String> cytoplasm = gProfilerResults.getProteinsWithGOTerm(GO_CYTOPLASM);
		final Set<String> organelle = gProfilerResults.getProteinsWithGOTerm(GO_ORGANELLE);

		// all together
		final VennData vennAll = getComplexesVennData("Cytoplasm vs Organelle (All complexes)", complexes, cytoplasm,
				organelle);
		final VennData vennKnown = getComplexesVennData("Cytoplasm vs Organelle (known complexes)", known, cytoplasm,
				organelle);
		final VennData vennUnknown = getComplexesVennData("Cytoplasm vs Organelle (unkown complexes)", unknown,
				cytoplasm, organelle);
		final File fig1File = File.createTempFile("fig1B", ".tmp");
		fig1File.deleteOnExit();
		final FileWriter fw = new FileWriter(fig1File);
		fw.write(vennAll.toString() + "\n\n\n");
		fw.write(vennKnown.toString() + "\n\n\n");
		fw.write(vennUnknown.toString());
		fw.close();
		FileUtils.separatedValuesToXLSX(fig1File.getAbsolutePath(), getExcelFile().getAbsolutePath(), "\t", "Fig1B");

	}

	private VennData getComplexesVennData(String vennDataTitle, List<ProteinComplex> complexes, Set<String> cytoplasm,
			Set<String> organelle) throws IOException {
		final List<ProteinComplex> cyto = new ArrayList<ProteinComplex>();
		final List<ProteinComplex> orga = new ArrayList<ProteinComplex>();
		for (final ProteinComplex proteinComplex : complexes) {
			int numCyto = 0;
			int numOrga = 0;
			for (final ProteinComponent component : proteinComplex.getComponentList()) {
				if (cytoplasm.contains(component.getAcc())) {
					numCyto++;
				}
				if (organelle.contains(component.getAcc())) {
					numOrga++;
				}
			}
			if (numCyto > numOrga) {
				cyto.add(proteinComplex);
			} else if (numOrga > numCyto) {
				orga.add(proteinComplex);
			} else {
				cyto.add(proteinComplex);
				orga.add(proteinComplex);
			}
		}
		final VennData ret = new VennData(vennDataTitle, "Cytoplasm", cyto, "Organelle", orga, null, null);
		return ret;
	}

	private File getExcelFile() {
		return new File(this.epicFolder.getAbsolutePath() + File.separator + getEpicFolderName() + "_figures.xlsx");
	}

	private File getEpicTotalProteinsGOCCFile() {
		return new File(
				this.epicFolder.getAbsolutePath() + File.separator + getEpicFolderName() + "_Total_Proteins_GO_CC.csv");
	}

	private File getEpicPPIsGOCCFile() {
		return new File(this.epicFolder.getAbsolutePath() + File.separator + getEpicFolderName() + "_ppi_GO_CC.csv");
	}

	private File getEpicPPIsFile() {
		return new File(this.epicFolder.getAbsolutePath() + File.separator + "Out.rf.exp.pred.txt");
	}

	private String getEpicFolderName() {
		return FilenameUtils.getBaseName(this.epicFolder.getAbsolutePath());
	}

	public FiguresGenerator(File epicFolder, String experimentFile) {
		this.projectSummaryFile = new File(experimentFile);
		this.experimentName = FilenameUtils.getBaseName(projectSummaryFile.getAbsolutePath());
		this.epicFolder = epicFolder;
	}

	public SeparationExperiment getSeparationExperiment() throws IOException {
		if (separationExperiment == null) {
			this.separationExperiment = ProteinComplexAnalyzer.loadProjectSummaryFileNEW(experimentName,
					projectSummaryFile);
		}
		return this.separationExperiment;
	}

}
