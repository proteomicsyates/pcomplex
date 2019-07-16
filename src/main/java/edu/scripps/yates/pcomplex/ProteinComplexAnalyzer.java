package edu.scripps.yates.pcomplex;

import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.proteored.pacom.analysis.charts.HeatChart;
import org.proteored.pacom.analysis.charts.HeatMapChart;

import com.ipa.ip2.api.exception.APIException;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.SftpException;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.pcomplex.db.ComplexPortalDB;
import edu.scripps.yates.pcomplex.db.CoreCorumDB;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.ftp.MySftpProgressMonitor;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComplexExistenceCriteria;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import edu.scripps.yates.pcomplex.util.PComplexUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.sequence.PeptideSequenceProperties;
import edu.scripps.yates.utilities.util.Pair;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinComplexAnalyzer {
	private final static Logger log = Logger.getLogger(ProteinComplexAnalyzer.class);
	public static final String[] projectsNames = { //
			//// "22_08_18_HEK_293_PCP" //
			//// , "22_10_12_HEK_293_PCP"//
			//// , "EvoSep"//
			// // ,
			//// , "150mM NaCl Size Exclusion PCP"//
			//// , "WCX_1" //
			//// ,
			"2019_05_17_SEC", //
			// "PCP_EvoSep_Elite_SEC1", // **
			// "PCP_EvoSep_VelosPro_WCX1", // **
			//
			// "WCX_SEC_Tandem", // **
			// "dSEC", // **
			// "2019_05_17_dSEC" **
	};

	// public static final String[] projectsNames = { //
	//
	// "FF_Drug_HEK" //

	// };
	private static boolean runOverlapping = false;

	public static String basePath = "Z:\\share\\Salva\\data\\asaviola\\protein complexes";
	public static String downloadFilesPath = basePath + "\\experiments";
	public static final String uniprotReleasesFolder = "Z:\\share\\Salva\\data\\uniprotKB";
	private static int defaultMinNumComponentsInComplex = 2;
	private final static DateFormat datef = new SimpleDateFormat("MM-dd-yyyy HH_mm_ss");
	public static String TAXONOMY;
	private static DecimalFormat df = new DecimalFormat("#.#");
	private final int minNumComponentsInComplex;
	private ProteinComplexExistenceCriteria existenceCriteria;
	private static final boolean downloadFiles = false;
	public static boolean useHUMAP = true;
	public static boolean useComplexPortalDB = false;
	public static boolean useCoreCorumDB = true;
	private static List<ProteinComplexDB> dbs;
	private final boolean getUniqueProteinsPerFraction = true;

	//
	private final boolean generateHeatMapsForAllComplexes = true;
	private final boolean compareWithProteinComplexDBs = true;
	private final boolean generateHeatMapsForIndividualProteinComplexes = false;

	private UniprotGeneMapping geneMapping;
	private FileWriter fileWriterForR;
	private boolean heatmapFolderCreated;

	private static Map<String, String> namesByProjectName = new THashMap<String, String>();
	private static String experimentNamePattern;

	public static void main(String[] args) {
		try {

			namesByProjectName.put("PCP_EvoSep_Elite_SEC1", "SEC");
			namesByProjectName.put("PCP_EvoSep_VelosPro_WCX1", "WCX");
			namesByProjectName.put("WCX_SEC_Tandem", "WCX_SEC");
			// basePath = "Z:\\share\\Salva\\data\\jim\\proteinComplexes";
			downloadFilesPath = basePath + "\\experiments";
			experimentNamePattern = "Cond_D_Rep_4_";
			final ProteinComplexAnalyzer pcan = new ProteinComplexAnalyzer(defaultMinNumComponentsInComplex);
			// USING CRITERIA
			pcan.setProteinComplexExistenceCriteria(new ProteinComplexExistenceCriteria());
			final List<SeparationExperiment> exps = new ArrayList<SeparationExperiment>();
			for (final String projectsName : projectsNames) {
				final SeparationExperiment exp = pcan.run(new File(args[0]), projectsName, downloadFiles);
				exps.add(exp);
				for (final DataType dataType : DataType.values()) {
					final File out = exp.exportToCSV(new File(args[0]).getParentFile(), dataType);
					log.info("File created at: " + out.getAbsolutePath());
				}
			}
			pcan.getFileWriterForR().close();
			if (runOverlapping) {
				SeparationExperiment exp1 = null;
				SeparationExperiment exp2 = null;
				SeparationExperiment exp3 = null;
				if (exps.size() > 0) {
					exp1 = exps.get(0);
				}
				if (exps.size() > 1) {
					exp2 = exps.get(1);
				}
				if (exps.size() > 2) {
					exp3 = exps.get(2);
				}
				pcan.compareSeparationExperiments(exp1, exp2, exp3);
			}

		} catch (final Exception e) {
			e.printStackTrace();
			log.error("Some error occurred...exiting now.");
			System.exit(-1);
		}
		log.info("Everything ok. Program finished.");
		System.exit(0);
	}

	private void setProteinComplexExistenceCriteria(ProteinComplexExistenceCriteria proteinComplexExistenceCriteria) {
		existenceCriteria = proteinComplexExistenceCriteria;
	}

	public static List<ProteinComplexDB> getDBs() throws IOException {
		if (dbs == null) {
			// add to list of DBs
			dbs = new ArrayList<ProteinComplexDB>();
			if (useHUMAP) {
				final ProteinComplexDB proteinComplexesORG = loadHUMapProteinComplexes();
				log.info(proteinComplexesORG.getProteinComplexes().size() + " protein complexes in "
						+ proteinComplexesORG.getName());
				dbs.add(proteinComplexesORG);
			}
			if (useComplexPortalDB) {
				final ProteinComplexDB complexPortalDB = loadComplexPortalProteinComplexes();
				log.info(complexPortalDB.getProteinComplexes().size() + " protein complexes in "
						+ complexPortalDB.getName());
				dbs.add(complexPortalDB);
			}
			if (useCoreCorumDB) {
				final ProteinComplexDB coreCorumDB = loadCoreCorumProteinComplexes();
				log.info(coreCorumDB.getProteinComplexes().size() + " protein complexes in " + coreCorumDB.getName());
				dbs.add(coreCorumDB);
			}
		}
		return dbs;
	}

	private void compareSeparationExperiments(SeparationExperiment exp1, SeparationExperiment exp2,
			SeparationExperiment exp3) throws IOException {

		for (final ProteinComplexDB proteinComplexDB : getDBs()) {
			// get all complete protein complexes
			Set<ProteinComplex> completedProteinComplexes1 = null;
			if (existenceCriteria == null) {
				completedProteinComplexes1 = exp1.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex);
			} else {
				completedProteinComplexes1 = exp1.getCompleteComplexes(proteinComplexDB, existenceCriteria,
						minNumComponentsInComplex);
			}
			Set<ProteinComplex> completedProteinComplexes2 = null;
			if (existenceCriteria == null) {
				completedProteinComplexes2 = exp2.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex);
			} else {
				completedProteinComplexes2 = exp2.getCompleteComplexes(proteinComplexDB, existenceCriteria,
						minNumComponentsInComplex);
			}
			String exp3Name = null;
			Set<ProteinComplex> completedProteinComplexes3 = null;
			if (exp3 != null) {
				if (existenceCriteria == null) {
					completedProteinComplexes3 = exp3.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex);
				} else {
					completedProteinComplexes3 = exp3.getCompleteComplexes(proteinComplexDB, existenceCriteria,
							minNumComponentsInComplex);
				}
				exp3Name = exp3.getProjectName();
			}
			final String title = proteinComplexDB.getName();
			final VennData vennData = new VennData(title, exp1.getProjectName(), completedProteinComplexes1,
					exp2.getProjectName(), completedProteinComplexes2, exp3Name, completedProteinComplexes3);

			System.out.println(vennData);
		}
	}

	public ProteinComplexAnalyzer(int minNumComponentsInComplex2) {
		TAXONOMY = "HUMAN";
		minNumComponentsInComplex = minNumComponentsInComplex2;
	}

	private static String dateString;

	private SeparationExperiment run(File propertiesFile, String projectName, boolean downloadFiles)
			throws IOException, JSchException, SftpException, APIException {
		geneMapping = UniprotGeneMapping.getInstance(new File(uniprotReleasesFolder), TAXONOMY);

		final IP2Util ip2Util = new IP2Util(new MySftpProgressMonitor(System.out), propertiesFile, projectName);
		if (downloadFiles) {
			downloadFiles(ip2Util, projectName);
		}
		return analyzeFractions(projectName);
	}

	private SeparationExperiment analyzeFractions(String projectName) throws IOException {
		final File projectFolder = getProjectFolder(projectName);
		final File projectSummaryFile = getProjectSummaryFile(projectName);
		if (!projectSummaryFile.exists() || projectSummaryFile.length() == 0) {
			FileWriter fw = null;
			try {
				fw = new FileWriter(projectSummaryFile);
				final File[] dtaSelectFiles = projectFolder.listFiles();
				downloadUniprotAnnotations(dtaSelectFiles);
				if (dtaSelectFiles != null) {
					for (final File file : dtaSelectFiles) {
						if (!file.isFile()) {
							continue;
						}
						final String fileName = FilenameUtils.getBaseName(file.getAbsolutePath());
						final int fractionNumber = Integer.valueOf(fileName.substring(0, fileName.indexOf("-")));
						String fractionName = fileName.substring(fileName.indexOf("-") + 1);
						fractionName = fractionName.substring(0, fractionName.lastIndexOf("-"));
						fw.write(fractionNumber + "\t" + fractionName);
						final DTASelectParser parser = new DTASelectParser(file);
						parser.enableProteinMergingBySecondaryAccessions(getUPLR(), null);
						parser.setOnlyReadProteins(true);
						parser.setDecoyPattern("Reverse");
						final Set<String> proteinAccs = parser.getProteinMap().keySet();
						for (final String proteinACC : proteinAccs) {
							String geneName = "";
							final Map<String, Entry> annotatedProtein = getUPLR().getAnnotatedProtein(null, proteinACC);
							final StringBuilder others = new StringBuilder();
							Double mw = null;
							if (annotatedProtein != null && annotatedProtein.containsKey(proteinACC)) {
								final Entry entry = annotatedProtein.get(proteinACC);
								mw = UniprotEntryUtil.getMolecularWeightInDalton(entry);
								final List<Pair<String, String>> geneNames = UniprotEntryUtil.getGeneName(entry, false,
										true);
								if (!geneNames.isEmpty()) {
									for (int i = 0; i < geneNames.size(); i++) {
										if (i == 0) {
											geneName = geneNames.get(i).getFirstelement();
										} else {
											if (!"".equals(others.toString())) {
												others.append(" | ");
											}
											others.append(geneNames.get(i).getFirstelement());
										}
									}
								}
								final Set<String> ensgids = UniprotEntryUtil.getENSGIDs(entry);
								for (final String ensgid : ensgids) {
									if (!"".equals(others.toString())) {
										others.append(" | ");
									}
									others.append(ensgid);
								}
							}
							final Integer spc = parser.getProteinMap().get(proteinACC).getSpectrumCount();
							final float nsaf = parser.getProteinMap().get(proteinACC).getNsaf();
							fw.write("\t" + proteinACC + " | " + geneName + " | " + mw + " | " + others.toString()
									+ " | " + spc + " | " + nsaf);
						}
						fw.write("\n");
					}
				}
			} finally {
				if (fw != null) {
					fw.close();
				}
			}
		}
		final SeparationExperiment experiment = loadProjectSummaryFile(projectName, projectSummaryFile);
		log.info(experiment.getFractions().size() + " fractions in project " + experiment.getProjectName());

		// load DBs
		// add to list of DBs
		final List<ProteinComplexDB> dbs = getDBs();

		FileWriter fw = null;

		// writting file for R
		// header
		getFileWriterForR().write(
				"Experiment\tDB\tType\tComplexID\tComplexName\tProtein\tFraction\tMW_TH\tMW_EXP\tPI_TH\tPI_EXP\n");
		try {
			fw = new FileWriter(getComparisonReportFile(projectName), false);

			fw.write("Analyzing " + experiment.getFractions().size() + " fractions from experiment "
					+ experiment.getProjectName() + "\n");

			fw.write("Fraction #\t");
			final List<Fraction> sortedFractions = experiment.getSortedFractions();
			for (final Fraction fraction : sortedFractions) {
				fw.write(fraction.getFractionNumber() + "\t");
			}
			fw.write("\n");
			// keep sets of accessions of each fraction
			final Map<Fraction, Set<String>> proteinsPerFraction = new THashMap<Fraction, Set<String>>();

			// num proteins in fraction
			fw.write("Num proteins in fraction\t");
			final Set<Protein> totalProteins = new THashSet<Protein>();
			for (final Fraction fraction : sortedFractions) {
				fw.write(fraction.getProteins().size() + "\t");
				totalProteins.addAll(fraction.getProteins());
				//
				if (getUniqueProteinsPerFraction) {
					final Set<String> accs = fraction.getProteins().stream().map(p -> p.getAcc())
							.collect(Collectors.toSet());
					proteinsPerFraction.put(fraction, accs);
				}
			}

			if (getUniqueProteinsPerFraction) {
				// num proteins unique to each fraction
				fw.write("\nUNIQUE proteins in fraction\t");
				for (int fractionIndex = 0; fractionIndex < sortedFractions.size(); fractionIndex++) {
					final Fraction fraction = sortedFractions.get(fractionIndex);
					final Set<String> set = proteinsPerFraction.get(fraction);
					final List<String> proteins = new ArrayList<String>();
					proteins.addAll(set);
					for (int fractionIndex2 = 0; fractionIndex2 < sortedFractions.size(); fractionIndex2++) {
						if (fractionIndex2 != fractionIndex) {
							final Fraction fraction2 = sortedFractions.get(fractionIndex2);
							final Set<String> set2 = proteinsPerFraction.get(fraction2);
							final Iterator<String> iterator = proteins.iterator();
							while (iterator.hasNext()) {
								final String protein = iterator.next();
								if (set2.contains(protein)) {
									iterator.remove();
								}
							}
						}
					}
					fw.write(proteins.size() + "\t");
				}
			}
			fw.write("\nNumber of total proteins in all fractions:\t" + totalProteins.size() + "\n");
			// annotate all
			log.info("Getting annotations of experimental proteins from UniprotKB");
			getUPLR().getAnnotatedProteins(null,
					totalProteins.stream().map(p -> p.getAcc()).collect(Collectors.toSet()));
			//

			// number of proteins with transmembrane regions
			int totalNumTransmembraneProteins = 0;
			for (final Protein protein : totalProteins) {
				if (PComplexUtil.hasTransmembraneRegion(protein.getAcc(), getUPLR())) {
					totalNumTransmembraneProteins++;
				}
			}
			fw.write("Total number of proteins with transmembrane regions (from UniprotKB annotation):\t"
					+ totalNumTransmembraneProteins + "\n");
			final double percentage = 100.0 * totalNumTransmembraneProteins / totalProteins.size();
			fw.write("% proteins with transmembrane regions (among the total detected):\t" + df.format(percentage)
					+ " %\n");

			if (compareWithProteinComplexDBs) {
				for (final ProteinComplexDB proteinComplexDB : dbs) {
					final String message = "Comparing with " + proteinComplexDB.getName() + ", that contains "
							+ proteinComplexDB.getProteinComplexes().size() + " protein complexes:";
					log.info(message);
					fw.write("\n" + message + "\n");
					compareWithDB(experiment, proteinComplexDB, fw);
					fw.write("\n\n");
				}
			}

			// HEATMAPS

			for (final ProteinComplexDB proteinComplexDB : getDBs()) {
				if (generateHeatMapsForIndividualProteinComplexes) {
					log.info("Generating heat maps from " + proteinComplexDB.getName() + " ...");
					for (int numFraction = 0; numFraction < sortedFractions.size(); numFraction++) {
						final Set<ProteinComplex> complexesViewed = new THashSet<ProteinComplex>();
						final Fraction fraction = sortedFractions.get(numFraction);
						final Set<ProteinComplex> completeComplexesOfFraction = fraction
								.getCompleteComplexes(proteinComplexDB, existenceCriteria, 4);
						final ProgressCounter counter = new ProgressCounter(completeComplexesOfFraction.size(),
								ProgressPrintingType.PERCENTAGE_STEPS, 0);

						for (final ProteinComplex pComplex : completeComplexesOfFraction) {
							if (complexesViewed.contains(pComplex)) { // to
																		// avoid
																		// create
																		// the
																		// heatmap
																		// for
																		// the
																		// same
																		// complex,
																		// several
																		// times).
								continue;
							}
							complexesViewed.add(pComplex);
							counter.increment();
							final String print = counter.printIfNecessary();
							if (print != null && !"".equals(print)) {
								log.info(proteinComplexDB.getName() + "\tFraction " + fraction.getFractionNumber()
										+ " (" + (numFraction + 1) + "/" + sortedFractions.size() + ")" + "\t" + print);
							}
							generateProteinComplexHeatMap(pComplex, proteinComplexDB, experiment, QuantType.SPC);
						}
					}
				}
				if (generateHeatMapsForAllComplexes) {
					generateHeatMapOfAllProteinComplexes(proteinComplexDB, experiment, QuantType.SPC,
							HeatChart.SCALE_EXPONENTIAL);
					generateHeatMapOfAllProteinComplexes(proteinComplexDB, experiment, QuantType.SPC,
							HeatChart.SCALE_LINEAR);
					generateHeatMapOfAllProteinComplexes(proteinComplexDB, experiment, QuantType.SPC,
							HeatChart.SCALE_LOGARITHMIC);
				}
			}

		} finally {
			if (fw != null) {
				fw.close();
			}
		}

		return experiment;
	}

	private void downloadUniprotAnnotations(File[] dtaSelectFiles) throws IOException {
		final Set<String> proteinAccs = new THashSet<String>();

		if (dtaSelectFiles != null) {
			for (final File file : dtaSelectFiles) {
				if (!file.isFile()) {
					continue;
				}

				final DTASelectParser parser = new DTASelectParser(file);
				parser.enableProteinMergingBySecondaryAccessions(null, null);
				parser.setOnlyReadProteins(true);
				parser.setDecoyPattern("Reverse");
				proteinAccs.addAll(parser.getProteinMap().keySet());
			}
		}

		getUPLR().getAnnotatedProteins(null, proteinAccs);
	}

	/**
	 * This heatmap will have a proteinComplex in each row and the quantiattive
	 * value for the complex in each fraction in each column. The rows are
	 * sorted by that value
	 * 
	 * @param proteinComplexDB
	 * @param experiment
	 * @param quantType
	 * @param colorScale
	 *            from HeatChart.SCALE_LINEAR
	 * @throws IOException
	 */
	private void generateHeatMapOfAllProteinComplexes(ProteinComplexDB proteinComplexDB,
			SeparationExperiment experiment, QuantType quantType, double colorScale) throws IOException {
		final List<Fraction> sortedFractions = experiment.getSortedFractions();
		final List<ProteinComplex> complexes = new ArrayList<ProteinComplex>();
		if (existenceCriteria != null) {
			complexes.addAll(
					experiment.getCompleteComplexes(proteinComplexDB, existenceCriteria, minNumComponentsInComplex));
		} else {
			complexes.addAll(experiment.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex));
		}
		// sort by quant type (MW, SPC or NSAF)
		Collections.sort(complexes, new Comparator<ProteinComplex>() {

			@Override
			public int compare(ProteinComplex o1, ProteinComplex o2) {
				Double v1 = null;
				Double v2 = null;

				switch (quantType) {
				case MW:
					try {
						v1 = o1.getMW(geneMapping, getUPLR());
					} catch (final IOException e) {
						e.printStackTrace();
					}
					try {
						v2 = o2.getMW(geneMapping, getUPLR());
					} catch (final IOException e) {
						e.printStackTrace();
					}
					break;
				// case SPC:
				// try {
				// v1 = o1.getSPC(geneMapping, getUPLR());
				// } catch (final IOException e) {
				// e.printStackTrace();
				// }
				// try {
				// v2 = o2.getSPC(geneMapping, getUPLR());
				// } catch (final IOException e) {
				// e.printStackTrace();
				// }
				// break;
				// case NSAF:
				// try {
				// v1 = o1.getNSAF(geneMapping, getUPLR());
				// } catch (final IOException e) {
				// e.printStackTrace();
				// }
				// try {
				// v2 = o2.getNSAF(geneMapping, getUPLR());
				// } catch (final IOException e) {
				// e.printStackTrace();
				// }
				// break;
				default:
					break;
				}
				if (v1 != v2) {
					return Double.compare(v1, v2);
				} else {
					if (v1 == null) {
						return 1;
					} else {
						return -1;
					}
				}
			}
		});
		final String projectName = experiment.getProjectName();
		final double[][] dataset = new double[complexes.size()][sortedFractions.size()];
		final List<String> rowList = complexes.stream().map(complex -> complex.getId()).collect(Collectors.toList());
		final List<String> columnList = new ArrayList<String>();
		for (int numFraction = 0; numFraction < sortedFractions.size(); numFraction++) {
			final Fraction fraction = sortedFractions.get(numFraction);
			columnList.add(String.valueOf(numFraction + 1));
			Set<ProteinComplex> complexesInFraction = null;
			if (existenceCriteria != null) {
				complexesInFraction = fraction.getCompleteComplexes(proteinComplexDB, existenceCriteria,
						minNumComponentsInComplex);
			} else {
				complexesInFraction = fraction.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex);
			}
			for (final ProteinComplex proteinComplexInFraction : complexesInFraction) {
				final double mw = proteinComplexInFraction.getMW(geneMapping, getUPLR());
				final int index = complexes.indexOf(proteinComplexInFraction);
				dataset[index][numFraction] = mw;
			}
		}

		final String title = "MWs of detected complexes (" + proteinComplexDB.getName() + ")";
		final HeatMapChart chart = new HeatMapChart(title, dataset, rowList, columnList, colorScale);
		chart.setTitleFon(new Font("Verdana", Font.PLAIN, 5));
		final File imageFile = getHeatMapImageFileForAllComplexes(projectName, experiment.getProjectName(),
				proteinComplexDB.getName(), colorScale, quantType);
		chart.saveImage(imageFile);

	}

	/**
	 * The heatmap, per protein complex, will have in each row the component of
	 * the complex and in each column the quantitative value of the protein in a
	 * particular fraction
	 * 
	 * @param pComplex
	 * @param proteinComplexDB
	 * @param experiment
	 * @param numFraction
	 * @throws IOException
	 */
	private void generateProteinComplexHeatMap(ProteinComplex pComplex, ProteinComplexDB proteinComplexDB,
			SeparationExperiment experiment, QuantType quantType) throws IOException {
		final List<Fraction> sortedFractions = experiment.getSortedFractions();
		final String projectName = experiment.getProjectName();
		final List<String> components = pComplex.getComponentList();
		final double[][] dataset = new double[pComplex.getComponents().size()][sortedFractions.size()];
		final List<String> rowList = components;
		final List<String> columnList = new ArrayList<String>();
		for (int numFraction2 = 0; numFraction2 < sortedFractions.size(); numFraction2++) {
			final Fraction fraction2 = sortedFractions.get(numFraction2);
			columnList.add(String.valueOf(numFraction2 + 1));
			for (int numComponent = 0; numComponent < components.size(); numComponent++) {
				final String component = components.get(numComponent);
				final Set<Protein> proteinsInFraction = fraction2.getProteinByKey(component);
				double quantValue = 0;
				if (proteinsInFraction != null) {
					for (final Protein protein : proteinsInFraction) {
						if (quantType == QuantType.SPC) {
							quantValue += protein.getSpc();
						}
					}

				}
				dataset[numComponent][numFraction2] = quantValue;
			}
		}

		String title = "Complex: " + pComplex.getId();
		if (pComplex.getName() != null) {
			title += ":" + pComplex.getName();
		}
		title += " (" + components.size() + " components)";
		title += "  (" + proteinComplexDB.getName() + ")";
		final HeatMapChart chart = new HeatMapChart(title, dataset, rowList, columnList, HeatChart.SCALE_LINEAR);
		chart.setTitleFon(new Font("Verdana", Font.PLAIN, 5));
		final File imageFile = getHeatMapImageFileForProteinComplex(projectName, experimentNamePattern,
				proteinComplexDB.getName(), pComplex.getId(), pComplex.getComponents().size());
		chart.saveImage(imageFile);
	}

	private File getComparisonReportFile(String projectName) {
		final String pathname = basePath + File.separator + getDateString();
		final File parentFolder = new File(pathname);
		if (!parentFolder.exists()) {
			parentFolder.mkdirs();
		}
		return new File(pathname + File.separator + projectName + "_" + getDateString() + "_comparison_output.tsv");
	}

	private File getComparisonReportFileforR() {
		final String pathname = basePath + File.separator + getDateString();
		final File parentFolder = new File(pathname);
		if (!parentFolder.exists()) {
			parentFolder.mkdirs();
		}
		return new File(pathname + File.separator + "data_for_R.tsv");
	}

	private String getNameFromProjectName(String projectName) {
		if (namesByProjectName.containsKey(projectName)) {
			return namesByProjectName.get(projectName);
		}
		return projectName;
	}

	private File getHeatMapImageFileForProteinComplex(String projectName, String experimentName, String dbName,
			String pComplex, int numComponents) {
		final String pathname = basePath + File.separator + "heatmaps";
		final File parentFolder = new File(pathname);
		if (!heatmapFolderCreated && !parentFolder.exists()) {
			parentFolder.mkdirs();
			heatmapFolderCreated = true;
		}
		return new File(pathname + File.separator + projectName + "_" + experimentName + "_" + numComponents + "_"
				+ dbName + "_" + pComplex + ".png");
	}

	private File getHeatMapImageFileForAllComplexes(String projectName, String experimentName, String dbName,
			double colorScale, QuantType quantType) {
		final String pathname = basePath + File.separator + "heatmaps";
		final File parentFolder = new File(pathname);
		if (!parentFolder.exists()) {
			parentFolder.mkdirs();
		}
		String scale = "";
		if (colorScale == HeatChart.SCALE_EXPONENTIAL) {
			scale = "expScale";
		} else if (colorScale == HeatChart.SCALE_LINEAR) {
			scale = "linearScale";
		} else if (colorScale == HeatChart.SCALE_LOGARITHMIC) {
			scale = "logScale";
		}
		return new File(pathname + File.separator + projectName + "_" + experimentName + "_" + dbName + "_color" + scale
				+ "_" + quantType + ".png");
	}

	private File getProteinUniquenessReportFile(String projectName) {
		final String pathname = basePath + File.separator + getDateString();
		final File parentFolder = new File(pathname);
		if (!parentFolder.exists()) {
			parentFolder.mkdirs();
		}
		return new File(
				pathname + File.separator + projectName + "_" + getDateString() + "_protein_uniqueness_output.tsv");
	}

	private void compareWithDB(SeparationExperiment experiment, ProteinComplexDB proteinComplexDB, FileWriter fw)
			throws IOException {
		final FileWriter fileWriterForR = getFileWriterForR();
		final Set<String> genesnotMapped = new THashSet<String>();

		final Set<String> uniprots = new THashSet<String>();
		for (final ProteinComplex complex : proteinComplexDB.getProteinComplexes()) {
			for (final String component : complex.getComponentList()) {
				final String uniProtACC = FastaParser.getUniProtACC(component);
				if (uniProtACC != null) {
					uniprots.add(uniProtACC);
				} else {
					final Set<String> mapped = geneMapping.mapGeneToUniprotACC(component);
					uniprots.addAll(mapped);
				}

			}
		}
		log.info("Getting annotations of " + uniprots.size() + " proteins in protein complex database  "
				+ proteinComplexDB.getName());

		getUPLR().getAnnotatedProteins(null, uniprots);

		int numFraction = 0;
		for (final Fraction fraction : experiment.getSortedFractions()) {
			numFraction++;
			Set<ProteinComplex> complexes = null;
			if (existenceCriteria == null) {
				complexes = fraction.getCompleteComplexes(proteinComplexDB, minNumComponentsInComplex);
			} else {
				complexes = fraction.getCompleteComplexes(proteinComplexDB, existenceCriteria,
						minNumComponentsInComplex);
			}
			for (final ProteinComplex complex : complexes) {
				final TDoubleArrayList allsubunitsMWs = new TDoubleArrayList();
				final TDoubleArrayList experimentaSubunitsMWs = new TDoubleArrayList();
				final TDoubleArrayList allsubunitsPIs = new TDoubleArrayList();
				final TDoubleArrayList experimentaSubunitsPIs = new TDoubleArrayList();

				final String complexID = complex.getId();
				final String complexName = complex.getName();
				for (final String componentKey : complex.getComponentList()) {
					boolean isDetectedExperimentally = false;
					if (!fraction.getProteinByKey(componentKey).isEmpty()) {
						isDetectedExperimentally = true;
					}
					Set<String> uniprotAccs = null;
					final String uniprotACC = FastaParser.getUniProtACC(componentKey);
					if (uniprotACC != null) {
						uniprotAccs = new THashSet<String>();
						uniprotAccs.add(uniprotACC);
					} else {
						final Set<String> mapGeneToUniprotACC = geneMapping.mapGeneToUniprotACC(componentKey);
						uniprotAccs = new THashSet<String>();

						for (final String acc : mapGeneToUniprotACC) {
							final Entry entry = getUPLR().getAnnotatedProtein(null, acc).get(acc);
							if (UniprotEntryUtil.isSwissProt(entry)) {
								uniprotAccs.add(acc);
							}
						}
						if (uniprotAccs.isEmpty()) {
							uniprotAccs.addAll(mapGeneToUniprotACC);
						}
					}
					if (uniprotAccs.isEmpty()) {
						genesnotMapped.add(componentKey);
					}
					final TDoubleArrayList ambiguousEntriesMWs = new TDoubleArrayList();
					final TFloatArrayList ambiguousEntriesPIs = new TFloatArrayList();
					for (final String acc : uniprotAccs) {
						final Entry entry = getUPLR().getAnnotatedProtein(null, acc).get(acc);
						if (entry != null) {
							final Double mw = UniprotEntryUtil.getMolecularWeightInDalton(entry);
							if (mw != null) {
								ambiguousEntriesMWs.add(mw);
							}
							final Float pi = PeptideSequenceProperties
									.calculatepI(UniprotEntryUtil.getProteinSequence(entry));
							if (pi != null) {
								ambiguousEntriesPIs.add(pi);
							}
						}
					}
					if (isDetectedExperimentally) {
						if (!ambiguousEntriesMWs.isEmpty()) {
							experimentaSubunitsMWs.add(Maths.mean(ambiguousEntriesMWs));
						}
						if (!ambiguousEntriesPIs.isEmpty()) {
							experimentaSubunitsPIs.add(Maths.mean(ambiguousEntriesPIs));
						}
					}
					if (!ambiguousEntriesMWs.isEmpty()) {
						allsubunitsMWs.add(Maths.mean(ambiguousEntriesMWs));
					}
					if (!ambiguousEntriesPIs.isEmpty()) {
						allsubunitsPIs.add(Maths.mean(ambiguousEntriesPIs));
					}
				}
				fileWriterForR.write(getNameFromProjectName(experiment.getProjectName()) + "\t");
				fileWriterForR.write(proteinComplexDB.getName() + "\t");
				fileWriterForR.write("COMPLEX\t");
				fileWriterForR.write(complexID + "\t" + complexName + "\t" + "-" + "\t" + numFraction + "\t");
				fileWriterForR.write(allsubunitsMWs.sum() + "\t" + experimentaSubunitsMWs.sum() + "\t");
				if (!allsubunitsPIs.isEmpty()) {
					fileWriterForR.write(Maths.mean(allsubunitsPIs) + "\t");
				} else {
					fileWriterForR.write("\t");
				}
				if (!experimentaSubunitsPIs.isEmpty()) {
					fileWriterForR.write(Maths.mean(experimentaSubunitsPIs) + "\n");
				} else {
					fileWriterForR.write("\n");
				}

			}
		}
		log.info(genesnotMapped.size() + " genes from complexes that are not mapped!");
		for (final String genenotMapped : genesnotMapped) {
			log.info(genenotMapped);
		}
		// Num proteins in each fraction
		final Set<Protein> totalProteins = new THashSet<Protein>();
		final Set<String> totalProteinKeys = new THashSet<String>();
		for (final Fraction fraction : experiment.getSortedFractions()) {
			fw.write(fraction.getProteins().size() + "\t");
			totalProteins.addAll(fraction.getProteins());
			fraction.getProteins().forEach((p) -> {
				totalProteinKeys.add(p.getAcc());
				totalProteinKeys.add(p.getGene());
			});
		}

		fw.write("\n");
		// num complete complexes
		final Set<ProteinComplex> totalCompleteProteinComplexes = new THashSet<ProteinComplex>();
		if (existenceCriteria == null) {
			fw.write("Num complete complexes\t");
			for (final Fraction fraction : experiment.getSortedFractions()) {
				final Set<ProteinComplex> completeComplexes = proteinComplexDB
						.getCompleteProteinComplexes(fraction.getProteinsAndGenes(), minNumComponentsInComplex);
				totalCompleteProteinComplexes.addAll(completeComplexes);
				fw.write(completeComplexes.size() + "\t");
			}

			fw.write("\n");
		} else {
			// num complete complexes with criteria
			fw.write("Num complete complexes\t");
			for (final Fraction fraction : experiment.getSortedFractions()) {
				final Set<ProteinComplex> completeComplexes = proteinComplexDB.getCompleteProteinComplexes(
						fraction.getProteinsAndGenes(), existenceCriteria, minNumComponentsInComplex);
				totalCompleteProteinComplexes.addAll(completeComplexes);
				fw.write(completeComplexes.size() + "\t");
			}
			fw.write("\n");
		}
		// num proteins in complexes
		fw.write("Num proteins in complete complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<Protein> proteins = null;
			if (existenceCriteria == null) {
				proteins = proteinComplexDB.getProteinsInCompleteProteinComplexes(fraction.getProteins(),
						minNumComponentsInComplex);
			} else {
				proteins = proteinComplexDB.getProteinsInCompleteProteinComplexes(fraction.getProteins(),
						existenceCriteria, minNumComponentsInComplex);
			}
			fw.write(proteins.size() + "\t");
		}
		fw.write("\n");
		// % proteins in complexes
		fw.write("% proteins in complete complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<Protein> proteins = null;
			if (existenceCriteria == null) {
				proteins = proteinComplexDB.getProteinsInCompleteProteinComplexes(fraction.getProteins(),
						minNumComponentsInComplex);
			} else {
				proteins = proteinComplexDB.getProteinsInCompleteProteinComplexes(fraction.getProteins(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final double percentage = 100.0 * proteins.size() / fraction.getProteins().size();
			fw.write(df.format(percentage) + "\t");
		}
		fw.write("\n");
		fw.write("AVG # of components in complete complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TIntArrayList numComponentts = new TIntArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				numComponentts.add(proteinComplex.getComponents().size());
			}
			if (numComponentts.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(df.format(Maths.mean(numComponentts)) + "\t");

			}
		}
		fw.write("\n");
		fw.write("AVG MW of complete complexes (proteins detected)\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TDoubleArrayList molWeights = new TDoubleArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				final TDoubleArrayList subunitsMWs = new TDoubleArrayList();
				for (final String componentKey : proteinComplex.getComponents()) {
					final TDoubleArrayList subunitMWs = new TDoubleArrayList();
					final Set<Protein> proteins = fraction.getProteinByKey(componentKey);

					for (final Protein protein : proteins) {
						if (protein.getMw() != null) {
							subunitMWs.add(protein.getMw());
						}
					}
					if (!subunitMWs.isEmpty()) {
						subunitsMWs.add(Maths.mean(subunitMWs));
					}
				}
				if (!subunitsMWs.isEmpty()) {
					molWeights.add(subunitsMWs.sum());
				}
			}
			if (molWeights.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(df.format(Maths.mean(molWeights)) + "\t");
			}
		}
		fw.write("\n");
		fw.write("AVG MW of complete complexes (detected+non-detected proteins)\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TDoubleArrayList molWeights = new TDoubleArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				final TDoubleArrayList subunitMWs = new TDoubleArrayList();
				for (final String componentKey : proteinComplex.getComponents()) {
					Set<String> uniprotAccs = null;
					final String uniprotACC = FastaParser.getUniProtACC(componentKey);
					if (uniprotACC != null) {
						uniprotAccs = new THashSet<String>();
						uniprotAccs.add(uniprotACC);
					} else {
						uniprotAccs = geneMapping.mapGeneToUniprotACC(componentKey);
					}
					for (final String acc : uniprotAccs) {
						final Entry entry = getUPLR().getAnnotatedProtein(null, acc).get(acc);
						if (entry != null) {
							final Double mw = UniprotEntryUtil.getMolecularWeightInDalton(entry);
							if (mw != null) {
								subunitMWs.add(mw);
								break;
							}
						}
					}
				}
				if (!subunitMWs.isEmpty()) {
					molWeights.add(subunitMWs.sum());
				}
			}
			if (molWeights.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(df.format(Maths.mean(molWeights)) + "\t");
			}
		}
		fw.write("\n");
		fw.write("AVG PI of complete complexes (detected proteins)\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TFloatArrayList pis = new TFloatArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				final TFloatArrayList subunitsPIs = new TFloatArrayList();
				for (final String componentKey : proteinComplex.getComponents()) {
					final TFloatArrayList subunitPIs = new TFloatArrayList();
					final Set<Protein> proteins = fraction.getProteinByKey(componentKey);

					for (final Protein protein : proteins) {
						final Entry entry = getUPLR().getAnnotatedProtein(null, protein.getAcc()).get(protein.getAcc());
						if (entry != null) {
							final float pi = PeptideSequenceProperties
									.calculatepI(UniprotEntryUtil.getProteinSequence(entry));
							subunitPIs.add(pi);
						}
					}
					if (!subunitPIs.isEmpty()) {
						subunitsPIs.add(Maths.mean(subunitPIs));
					}
				}
				if (!subunitsPIs.isEmpty()) {
					pis.add(Maths.mean(subunitsPIs));
				}
			}
			if (pis.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(df.format(Maths.mean(pis)) + "\t");
			}
		}
		fw.write("\n");
		fw.write("AVG PI of complete complexes (detected + non detected proteins)\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TFloatArrayList pis = new TFloatArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				final TFloatArrayList subunitsPIs = new TFloatArrayList();
				for (final String componentKey : proteinComplex.getComponents()) {
					Set<String> uniprotAccs = null;
					final String uniProtACC = FastaParser.getUniProtACC(componentKey);
					if (uniProtACC != null) {
						uniprotAccs = new THashSet<String>();
						uniprotAccs.add(uniProtACC);
					} else {
						uniprotAccs = geneMapping.mapGeneToUniprotACC(componentKey);
					}
					for (final String acc : uniprotAccs) {
						final Entry entry = getUPLR().getAnnotatedProtein(null, acc).get(acc);

						if (entry != null) {
							final float pi = PeptideSequenceProperties
									.calculatepI(UniprotEntryUtil.getProteinSequence(entry));
							subunitsPIs.add(pi);
							break;
						}
					}
				}
				if (!subunitsPIs.isEmpty()) {
					pis.add(Maths.mean(subunitsPIs));
				}
			}
			if (pis.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(df.format(Maths.mean(pis)) + "\t");
			}
		}
		fw.write("\n");
		fw.write("MAX # of components in complete complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TIntArrayList numComponentts = new TIntArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				numComponentts.add(proteinComplex.getComponents().size());
			}
			if (numComponentts.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(numComponentts.max() + "\t");
			}
		}
		fw.write("\n");
		fw.write("MIN # of components in complete complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			final TIntArrayList numComponentts = new TIntArrayList();
			for (final ProteinComplex proteinComplex : completeComplexes) {
				numComponentts.add(proteinComplex.getComponents().size());
			}
			if (numComponentts.isEmpty()) {
				fw.write(0 + "\t");
			} else {
				fw.write(numComponentts.min() + "\t");
			}
		}
		fw.write("\n");
		// num partial complexes
		final Set<ProteinComplex> totalPartialProteinComplexes = new THashSet<ProteinComplex>();
		fw.write("Num partial complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Set<ProteinComplex> partialComplexes = proteinComplexDB
					.getPartiallyCompletedProteinComplexes(fraction.getProteins());
			totalPartialProteinComplexes.addAll(partialComplexes);
			fw.write(partialComplexes.size() + "\t");
		}
		fw.write("\n");

		// num partial complexes
		fw.write("Num proteins in partial complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Set<Protein> proteins = proteinComplexDB
					.getProteinsInPartiallyCompletedProteinComplexes(fraction.getProteins());
			fw.write(proteins.size() + "\t");
		}
		fw.write("\n");

		// num partial complexes
		fw.write("% proteins in partial complexes\t");
		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Set<Protein> proteins = proteinComplexDB
					.getProteinsInPartiallyCompletedProteinComplexes(fraction.getProteins());
			final double percentage = 100.0 * proteins.size() / fraction.getProteins().size();
			fw.write(df.format(percentage) + "\t");
		}

		fw.write("\n");

		fw.write("Number of total completed protein complexes:\t" + totalCompleteProteinComplexes.size() + "\n");
		fw.write("Number of total partial protein complexes:\t" + totalPartialProteinComplexes.size() + "\n");

		// number of proteins in partially covered complexes that have
		// transmembrane regions
		final Set<String> totalProteinsNotCoveredInPartialComplexes = new THashSet<String>();
		for (final ProteinComplex complex : totalPartialProteinComplexes) {
			List<String> accs = complex.getComponentsNotInProteinsAndGenes(totalProteinKeys);
			if (complex.isGeneNames()) {
				accs = accs.stream().map(acc -> acc.toString() + "_" + TAXONOMY).collect(Collectors.toList());
			}
			totalProteinsNotCoveredInPartialComplexes.addAll(accs);
		}

		getUPLR().getAnnotatedProteins(null, totalProteinsNotCoveredInPartialComplexes);
		int totalTransmembraneProteinsNotCoveredInPartialComplexes = 0;
		for (final String acc : totalProteinsNotCoveredInPartialComplexes) {
			if (PComplexUtil.hasTransmembraneRegion(acc, getUPLR())) {
				totalTransmembraneProteinsNotCoveredInPartialComplexes++;
			}
		}
		fw.write("Total number of proteins not covered in partially covered complexes:\t"
				+ totalProteinsNotCoveredInPartialComplexes.size() + "\n");
		fw.write(
				"Total number of proteins with transmembrane regions that are among the ones that are not covered in partially covered complexes:\t"
						+ totalTransmembraneProteinsNotCoveredInPartialComplexes + "\n");

		final double percentage = 100.0 * totalTransmembraneProteinsNotCoveredInPartialComplexes
				/ totalProteinsNotCoveredInPartialComplexes.size();
		fw.write(
				"% of proteins with transmembrane regions that are among the ones that are not covered in partially covered complexes:\t"
						+ df.format(percentage) + " %\n");

		// write individual files per fraction
		writeIndividualFilesPerFraction(experiment, proteinComplexDB);

	}

	private FileWriter getFileWriterForR() throws IOException {
		if (fileWriterForR == null) {
			fileWriterForR = new FileWriter(getComparisonReportFileforR());
		}
		return fileWriterForR;
	}

	private void writeIndividualFilesPerFraction(SeparationExperiment experiment, ProteinComplexDB proteinComplexDB)
			throws IOException {
		for (final Fraction fraction : experiment.getSortedFractions()) {
			Set<ProteinComplex> completeComplexes = null;
			if (existenceCriteria == null) {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						minNumComponentsInComplex);
			} else {
				completeComplexes = proteinComplexDB.getCompleteProteinComplexes(fraction.getProteinsAndGenes(),
						existenceCriteria, minNumComponentsInComplex);
			}
			// print file with protein complexes for that fraction
			final FileWriter fwFraction = new FileWriter(getFractionCompleteProteinComplexesFile(
					experiment.getProjectName(), fraction, proteinComplexDB, false), false);
			fwFraction.write(completeComplexes.size() + " complexes are complete in fraction "
					+ fraction.getFractionNumber() + " " + fraction.getName() + "\n");
			for (final ProteinComplex proteinComplex : completeComplexes) {
				fwFraction.write(proteinComplex.toString() + "\n\n");
			}
			fwFraction.close();
		}

		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Set<ProteinComplex> partialComplexes = proteinComplexDB
					.getPartiallyCompletedProteinComplexes(fraction.getProteins());

			// print file with protein complexes for that fraction
			final FileWriter fwFraction = new FileWriter(getFractionPartialProteinComplexesFile(
					experiment.getProjectName(), fraction, proteinComplexDB, false), false);
			fwFraction.write(partialComplexes.size() + " complexes are partially present in fraction "
					+ fraction.getFractionNumber() + " " + fraction.getName() + "\n");
			for (final ProteinComplex complex : partialComplexes) {
				fwFraction.write(complex.getId() + " - " + complex.getName() + "\n");
				fwFraction.write("Present:\n");
				final List<String> presentComponents = complex
						.getComponentsInProteinsAndGenes(fraction.getProteinsAndGenes());
				for (final String prot : presentComponents) {
					fwFraction.write(prot.toString() + "+\t");

				}
				fwFraction.write("\nNot present:\n");
				final List<String> nonPresentComponents = complex
						.getComponentsNotInProteinsAndGenes(fraction.getProteinsAndGenes());
				for (final String prot : nonPresentComponents) {
					fwFraction.write(prot.toString() + "\t");
				}
				fwFraction.write("\n\n");
			}
			fwFraction.close();
		}
	}

	private File getFractionCompleteProteinComplexesFile(String experimentName, Fraction fraction,
			ProteinComplexDB proteinComplexDB, boolean b) {
		final File file = new File(basePath + File.separator + getDateString() + File.separator + experimentName
				+ File.separator + "Fr_" + "_" + fraction.getFractionNumber() + "_complete_complexes_"
				+ proteinComplexDB.getName() + ".tsv");
		if (!file.getParentFile().exists()) {
			file.getParentFile().mkdirs();
		}
		return file;
	}

	private File getFractionPartialProteinComplexesFile(String experimentName, Fraction fraction,
			ProteinComplexDB proteinComplexDB, boolean b) {
		final File file = new File(basePath + File.separator + getDateString() + File.separator + experimentName
				+ File.separator + "Fr_" + "_" + fraction.getFractionNumber() + "_partial_complexes_"
				+ proteinComplexDB.getName() + ".tsv");
		if (!file.getParentFile().exists()) {
			file.getParentFile().mkdirs();
		}
		return file;
	}

	public static ProteinComplexDB loadCoreCorumProteinComplexes() throws IOException {
		final ProteinComplexDB proteinComplexDB = new CoreCorumDB(getCoreCorumFile(), "Human");
		return proteinComplexDB;
	}

	private static File getCoreCorumFile() {
		final String path = basePath + File.separator + "corum" + File.separator + "coreComplexes.txt";
		return new File(path);
	}

	public static ProteinComplexDB loadComplexPortalProteinComplexes() throws IOException {
		final ProteinComplexDB proteinComplexDB = new ComplexPortalDB(getComplexPortalFile());
		return proteinComplexDB;
	}

	public static ProteinComplexDB loadHUMapProteinComplexes() throws IOException {
		final ProteinComplexDB proteinComplexDB = new ProteinComplexDB(getHUMapFile(), "Hu.MAP", getHUMapMappingFile(),
				true);
		return proteinComplexDB;
	}

	protected static SeparationExperiment loadProjectSummaryFile(String projectName, File projectSummaryFile)
			throws IOException {
		final SeparationExperiment ret = new SeparationExperiment(projectName);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(projectSummaryFile));
			String line = br.readLine();

			while (line != null) {
				final String[] split = line.split("\t");
				int fractionNumber = -1;
				String fractionName = null;
				for (int index = 0; index < split.length; index++) {
					final String tmp = split[index];
					if (index == 0) {
						fractionNumber = Integer.valueOf(tmp);
					} else if (index == 1) {
						fractionName = tmp;

					} else {
						final String[] split2 = tmp.split("\\|");
						final String proteinAcc = split2[0].trim();
						final String gene = split2[1].trim();
						Double mw = null;
						if (!"".equals(split2[2].trim()) && !"null".equals(split2[2].trim())) {
							mw = Double.valueOf(split2[2].trim());
						}
						final int spc = Integer.valueOf(split2[split2.length - 2].trim());
						final float nsaf = Float.valueOf(split2[split2.length - 1].trim());
						final Protein protein = new Protein(proteinAcc, gene, mw, spc, nsaf, fractionName);
						for (int i = 3; i < split2.length - 2; i++) {
							final String other = split2[i].trim();
							protein.addOther(other);
						}
						ret.addProtein(fractionName, fractionNumber, protein);
					}
				}
				line = br.readLine();
			}
			return ret;
		} finally {
			if (br != null) {
				br.close();
			}
		}
	}

	private UniprotProteinLocalRetriever getUPLR() {
		return new UniprotProteinLocalRetriever(new File(uniprotReleasesFolder), true);
	}

	private void downloadFiles(IP2Util ip2Util, String projectName) throws IOException, JSchException, SftpException {
		final Map<String, String> dtaSelectsInProject = ip2Util.getDTASelectsInProject();

		for (final String experimentName : dtaSelectsInProject.keySet()) {
			try {
				final int fraction = IP2Util.getFractionNumberFromExperimentName(experimentName, experimentNamePattern);

				if (fraction < 1) {
					continue;
				}
				log.info(fraction + "-\t" + experimentName);
				final File outputFile = getexperimentalDataOutputFile(projectName, experimentName,
						experimentNamePattern);
				final OutputStream outputStream = new FileOutputStream(outputFile, false);
				ip2Util.download(experimentName, dtaSelectsInProject.get(experimentName), outputStream);
				final DTASelectParser parser = new DTASelectParser(outputFile);
				final List<String> commandLineParameterStrings = parser.getCommandLineParameterStrings();
				for (final String param : commandLineParameterStrings) {
					if (param.trim().startsWith("-p")) {
						if (param.trim().startsWith("-p 2")) {
							// discarding that file
							outputFile.delete();
							log.info(outputFile.getAbsolutePath()
									+ " deleted because it has not correct search parameters (-p 2)");
						}
					}
				}

			} catch (final IllegalArgumentException e) {
				if (experimentNamePattern == null) {
					e.printStackTrace();
				}
			}
		}
	}

	private File getProjectFolder(String projectName) {
		final String folderPath = downloadFilesPath + File.separator + projectName;
		return new File(folderPath);
	}

	private File getProjectSummaryFile(String projectName) {
		final String path = downloadFilesPath + File.separator + projectName + ".tsv";
		return new File(path);
	}

	private static File getComplexPortalFile() {
		final String path = basePath + File.separator + "ComplexPortal" + File.separator + "homo_sapiens.tsv";
		return new File(path);
	}

	private static File getHUMapFile() {
		final String path = basePath + File.separator + "ProteinComplexes.org" + File.separator
				+ "genename_clusters.txt";
		return new File(path);
	}

	private static File getHUMapMappingFile() {
		final String path = basePath + File.separator + "ProteinComplexes.org" + File.separator + "nodeTable.txt";
		return new File(path);
	}

	private File getexperimentalDataOutputFile(String projectName, String experimentName,
			String experimentNamePattern) {
		final int fractionNumber = IP2Util.getFractionNumberFromExperimentName(experimentName, experimentNamePattern);
		final File projectFolder = getProjectFolder(projectName);
		projectFolder.mkdirs();
		final String path = projectFolder.getAbsolutePath() + File.separator + fractionNumber + "-" + experimentName
				+ "-DTASelect.txt";
		return new File(path);
	}

	private static String getDateString() {
		if (dateString == null) {
			dateString = datef.format(Calendar.getInstance().getTime());
		}
		return dateString;
	}
}
