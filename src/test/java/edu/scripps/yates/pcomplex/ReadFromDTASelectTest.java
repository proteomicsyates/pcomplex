package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.pcomplex.db.CoreCorumDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComplexExistenceCriteria;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class ReadFromDTASelectTest {

	private final File corumFile = new File(
			"Z:\\share\\Salva\\data\\asaviola\\protein complexes\\corum\\coreComplexes.txt");
	private final String[] corumOrganisms = null;// "Rat";
	private final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
	private final boolean filterCorumByPurificationMethod = false;
	private final int maxComponentsToBeSmall = ProteinComplexExistenceCriteria.default_MAX_COMPONENTS_TO_BE_SMALL;
	private final double percentageInBigComplexes = ProteinComplexExistenceCriteria.defaultPercentageInBigComplexes;
	private final double percentageInSmallComplexes = ProteinComplexExistenceCriteria.defaultPercentageInSmallComplexes;
	private final String folder = "Z:\\share\\Salva\\data\\Daniel\\ComplexAnalysis";
	private final String[] prefixes = { "ControlPCP", "Controlsal", "GAPpcp", "GAPsal", "GRM5pcp", "GRM5sal", "GSKpcp",
			"GSKsal", "MEKpcp", "MEKsal", "NR2Bpcp", "NRB2Bsal", "PPPpcp", "PPPsal", "STXpcp", "STXsal", "SYTpcp",
			"SYTsal" };
	private final String[] experimentKeys = { "gap", "grm5", "gsk", "mek", "nr2b", "ppp", "stx", "syt" };
	private final String[] salpcp = { "sal", "pcp" };
	private final int minNumComponents = 3;

	@Test
	public void readComplexesFromDTASelect() {
		final int numReps = 3;
		final List<String> fileNames = new ArrayList<String>();
		// controls
		for (final String salorpcp : salpcp) {
			for (final String experimentKey : experimentKeys) {

				fileNames.add("Control" + salorpcp + "%" + "_" + experimentKey);

			}
		}

		// the rest
		for (final String experimentKey : experimentKeys) {
			for (final String salorpcp : salpcp) {
				fileNames.add(experimentKey.toUpperCase() + salorpcp + "%");
			}
		}

		if (corumOrganisms != null && corumOrganisms.length > 0) {
			ProteinComplexAnalyzer.ORGANISMS = corumOrganisms; // take all and work with gene names

		}
		FileWriter fw = null;
		try {
			fw = new FileWriter(folder + File.separator + "complexes.txt");
			final CoreCorumDB corum = loadCorum();
			fw.write("exp\tNum complexes\n");

			for (final String fileName : fileNames) {
				final String fileName2 = folder + File.separator + fileName.replace("%", "") + "_complexes.txt";
				final FileWriter fw2 = new FileWriter(fileName2);
				final Set<String> proteinsAndGenes = new THashSet<String>();
				final TIntObjectHashMap<Set<String>> proteinsAndGenesByRep = new TIntObjectHashMap<Set<String>>();
				final TIntObjectHashMap<Set<String>> proteinComplexesIDsByRep = new TIntObjectHashMap<Set<String>>();
				for (int rep = 1; rep <= numReps; rep++) {
					final File dtaSelectFile = new File(
							folder + File.separator + fileName.replace("%", rep + "") + "_DTASelect-filter.txt");
					final List<ProteinComponent> components = new ArrayList<ProteinComponent>();
					DTASelectParser parser;

					parser = new DTASelectParser(dtaSelectFile);

					parser.setOnlyReadProteins(true);
					parser.setDecoyPattern("Reverse");
					for (final Protein protein : parser.getProteins()) {
						String gene = null;
						if (protein.getGenes() != null && !protein.getGenes().isEmpty()) {
							gene = protein.getGenes().iterator().next().getGeneID();
						}
						final ProteinComponent component = new ProteinComponent(protein.getAccession(), gene);
						components.add(component);
					}
					final Set<String> proteinsAndGenesrep = getProteinsAndGenes(components);
					proteinsAndGenesByRep.put(rep, proteinsAndGenesrep);
					proteinsAndGenes.addAll(proteinsAndGenesrep);
					final Set<ProteinComplex> proteinComplexes = corum.getCompleteProteinComplexes(proteinsAndGenesrep,
							getExistenceCriteria(), minNumComponents);
					proteinComplexesIDsByRep.put(rep,
							proteinComplexes.stream().map(c -> c.getId()).collect(Collectors.toSet()));

					fw.write(fileName.replace("%", rep + "") + "\t" + proteinComplexes.size() + "\n");
				}
				final Set<ProteinComplex> proteinComplexes = corum.getCompleteProteinComplexes(proteinsAndGenes,
						getExistenceCriteria(), minNumComponents);
				fw.write("Total in " + fileName.replace("%", "") + "\t" + proteinComplexes.size() + "\n");
				fw2.write("Total in " + fileName.replace("%", "") + "\t" + proteinComplexes.size() + "\n");
				fw2.write("ComplexID\tComplex name\t# componets in DB\tOrganism\t");
				for (int i = 1; i <= numReps; i++) {
					fw2.write("in rep " + i + "\t");
				}
				for (int i = 1; i <= numReps; i++) {
					fw2.write("# components in rep " + i + "\t");
				}
				fw2.write("\n");
				for (final ProteinComplex proteinComplex : proteinComplexes) {
					fw2.write(proteinComplex.getId() + "\t" + proteinComplex.getName() + "\t"
							+ proteinComplex.getComponentList().size() + "\t" + proteinComplex.getOrganism() + "\t");
					// is present or not
					for (int i = 1; i <= numReps; i++) {
						final Set<String> proteinComplexesIDs = proteinComplexesIDsByRep.get(i);
						if (!proteinComplexesIDs.contains(proteinComplex.getId())) {
							fw2.write("0");
						} else {
							fw2.write("1");
						}
						fw2.write("\t");
					}
					// # units
					for (int i = 1; i <= numReps; i++) {
						final Set<String> proteinsAndGenesRep = proteinsAndGenesByRep.get(i);
						final Set<String> proteinComplexesIDs = proteinComplexesIDsByRep.get(i);
						if (!proteinComplexesIDs.contains(proteinComplex.getId())) {
							fw2.write("-");
						} else {
							int num = 0;
							for (final ProteinComponent component : proteinComplex.getComponentList()) {
								if (proteinsAndGenesRep.contains(component.getAcc())
										|| proteinsAndGenesRep.contains(component.getGene())) {
									num++;
								}
							}
							fw2.write(String.valueOf(num));
						}
						fw2.write("\t");

					}
					fw2.write("\n");
				}
				fw2.close();
			}
			fw.write("\nCriteria are:\n" + getExistenceCriteria());
			fw.write("Protein complexes of " + (minNumComponents - 1) + " or less are discarded \n");
			fw.write("Comparing against CorumDB (" + corum.getProteinComplexes().size() + " complexes for "
					+ StringUtils.getSortedSeparatedValueStringFromChars(corumOrganisms, "-") + ")\n");
			fw.write("-----------------------------------------------------------------------------------\n");
			fw.close();
			saveToExcel(folder, "");
		} catch (final FileNotFoundException e) {
			e.printStackTrace();
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {

		}
	}

	@Test
	public void readComplexesFromDTASelect2() {
		final int numReps = 3;
		final String[] fileNames = { "PCP%_lysate", "Saline%_lysate" };

		if (corumOrganisms != null && corumOrganisms.length > 0) {
			ProteinComplexAnalyzer.ORGANISMS = corumOrganisms; // take all and work with gene names

		}
		FileWriter fw = null;
		try {
			fw = new FileWriter(folder + File.separator + "complexes.txt");
			final CoreCorumDB corum = loadCorum();
			fw.write("exp\tNum complexes\n");

			for (final String fileName : fileNames) {
				final String fileName2 = folder + File.separator + fileName.replace("%", "") + "_complexes.txt";
				final FileWriter fw2 = new FileWriter(fileName2);
				final Set<String> proteinsAndGenes = new THashSet<String>();
				final TIntObjectHashMap<Set<String>> proteinsAndGenesByRep = new TIntObjectHashMap<Set<String>>();
				final TIntObjectHashMap<Set<String>> proteinComplexesIDsByRep = new TIntObjectHashMap<Set<String>>();
				for (int rep = 1; rep <= numReps; rep++) {
					final File dtaSelectFile = new File(
							folder + File.separator + fileName.replace("%", rep + "") + "_DTASelect-filter.txt");
					final List<ProteinComponent> components = new ArrayList<ProteinComponent>();
					DTASelectParser parser;

					parser = new DTASelectParser(dtaSelectFile);

					parser.setOnlyReadProteins(true);
					parser.setDecoyPattern("Reverse");
					for (final Protein protein : parser.getProteins()) {
						String gene = null;
						if (protein.getGenes() != null && !protein.getGenes().isEmpty()) {
							gene = protein.getGenes().iterator().next().getGeneID();
						}
						final ProteinComponent component = new ProteinComponent(protein.getAccession(), gene);
						components.add(component);
					}
					final Set<String> proteinsAndGenesrep = getProteinsAndGenes(components);
					proteinsAndGenesByRep.put(rep, proteinsAndGenesrep);
					proteinsAndGenes.addAll(proteinsAndGenesrep);
					final Set<ProteinComplex> proteinComplexes = corum.getCompleteProteinComplexes(proteinsAndGenesrep,
							getExistenceCriteria(), minNumComponents);
					proteinComplexesIDsByRep.put(rep,
							proteinComplexes.stream().map(c -> c.getId()).collect(Collectors.toSet()));

					fw.write(fileName.replace("%", rep + "") + "\t" + proteinComplexes.size() + "\n");
				}
				final Set<ProteinComplex> proteinComplexes = corum.getCompleteProteinComplexes(proteinsAndGenes,
						getExistenceCriteria(), minNumComponents);
				fw.write("Total in " + fileName.replace("%", "") + "\t" + proteinComplexes.size() + "\n");
				fw2.write("Total in " + fileName.replace("%", "") + "\t" + proteinComplexes.size() + "\n");
				fw2.write("ComplexID\tComplex name\t# componets in DB\tOrganism\t");
				for (int i = 1; i <= numReps; i++) {
					fw2.write("in rep " + i + "\t");
				}
				for (int i = 1; i <= numReps; i++) {
					fw2.write("# components in rep " + i + "\t");
				}
				fw2.write("\n");
				for (final ProteinComplex proteinComplex : proteinComplexes) {
					fw2.write(proteinComplex.getId() + "\t" + proteinComplex.getName() + "\t"
							+ proteinComplex.getComponentList().size() + "\t" + proteinComplex.getOrganism() + "\t");
					// is present or not
					for (int i = 1; i <= numReps; i++) {
						final Set<String> proteinComplexesIDs = proteinComplexesIDsByRep.get(i);
						if (!proteinComplexesIDs.contains(proteinComplex.getId())) {
							fw2.write("0");
						} else {
							fw2.write("1");
						}
						fw2.write("\t");
					}
					// # units
					for (int i = 1; i <= numReps; i++) {
						final Set<String> proteinsAndGenesRep = proteinsAndGenesByRep.get(i);
						final Set<String> proteinComplexesIDs = proteinComplexesIDsByRep.get(i);
						if (!proteinComplexesIDs.contains(proteinComplex.getId())) {
							fw2.write("-");
						} else {
							int num = 0;
							for (final ProteinComponent component : proteinComplex.getComponentList()) {
								if (proteinsAndGenesRep.contains(component.getAcc())
										|| proteinsAndGenesRep.contains(component.getGene())) {
									num++;
								}
							}
							fw2.write(String.valueOf(num));
						}
						fw2.write("\t");

					}
					fw2.write("\n");
				}
				fw2.close();
			}
			fw.write("\nCriteria are:\n" + getExistenceCriteria());
			fw.write("Protein complexes of " + (minNumComponents - 1) + " or less are discarded \n");
			fw.write("Comparing against CorumDB (" + corum.getProteinComplexes().size() + " complexes for "
					+ StringUtils.getSortedSeparatedValueStringFromChars(corumOrganisms, "-") + ")\n");
			fw.write("-----------------------------------------------------------------------------------\n");
			fw.close();
			saveToExcel(folder, "PCP_Saline");
		} catch (final FileNotFoundException e) {
			e.printStackTrace();
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {

		}
	}

	private void saveToExcel(String folder, String fileName) throws IOException {
		final File[] listFiles = new File(folder).listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.contains("_complexes.txt")) {
					return true;
				}
				return false;
			}
		});
		final List<File> list = new ArrayList<File>();
		for (final File file : listFiles) {
			list.add(file);
		}
		Collections.sort(list, new Comparator<File>() {

			@Override
			public int compare(File o1, File o2) {
				return FilenameUtils.getBaseName(o1.getAbsolutePath())
						.compareTo(FilenameUtils.getBaseName(o2.getAbsolutePath()));
			}
		});
		final String outputXlsFilePath = folder + File.separator + fileName + "_complexes.xlsx";
		for (final File file : listFiles) {
			if (FilenameUtils.getBaseName(file.getAbsolutePath()).equals("complexes")) {
				System.out.println(file.getAbsolutePath());
			}
			final String sheetName = FilenameUtils.getBaseName(file.getAbsolutePath().replace("_complexes", ""));
			FileUtils.separatedValuesToXLSX(file.getAbsolutePath(), outputXlsFilePath, "\t", sheetName);
		}

	}

	private ProteinComplexExistenceCriteria getExistenceCriteria() {
		final ProteinComplexExistenceCriteria ret = new ProteinComplexExistenceCriteria();
		ret.setMaxComponentsToBeSmall(maxComponentsToBeSmall);
		ret.setPercentageInBigComplexes(percentageInBigComplexes);
		ret.setPercentageInSmallComplexes(percentageInSmallComplexes);
		return ret;
	}

	private Set<String> getProteinsAndGenes(List<ProteinComponent> components) {
		final Set<String> ret = new HashSet<String>();
		for (final ProteinComponent proteinComponent : components) {
			ret.add(proteinComponent.getAcc());
			if (proteinComponent.getGene() != null) {
				ret.add(proteinComponent.getGene());
			}
		}
		return ret;
	}

	private CoreCorumDB loadCorum() throws IOException {
		final CoreCorumDB corum = new CoreCorumDB(corumFile, corumOrganisms, getUPLR(),
				filterCorumByPurificationMethod);
		return corum;
	}

	private UniprotProteinLocalRetriever getUPLR() {
		return new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
	}
}
