package edu.scripps.yates.pcomplex.epic;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.ProteinProteinInteraction;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.set.hash.THashSet;

public class EpicResultsComparations {
	@Test
	public void compareToComplexDBs() {

		ProteinComplexAnalyzer.useComplexPortalDB = true;
		ProteinComplexAnalyzer.useCoreCorumDB = true;
		ProteinComplexAnalyzer.useHUMAP = true;
		final List<ProteinComplexDB> dBs = ProteinComplexAnalyzer.getDBs();
		final File folder = new File("C:\\Users\\salvador\\Desktop\\epic\\input");
		final double minOverlapScore = 0.25;
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			final File outputFile = epicComparator.compareWithDBs(dBs);
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());
			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Test
	public void compareToOtherEPICResult() {

		final File folder = new File("C:\\Users\\salvador\\Desktop\\epic\\input\\Kirkwood_et_al_out");
		final double minOverlapScore = 0.25;
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			final File outputFile = epicComparator.compareWithOtherEPICResult(
					new File("C:\\Users\\salvador\\Desktop\\epic\\input\\2019_11_19_Mixed_bed_SEC_SPC_out"),
					new File("C:\\Users\\salvador\\Desktop\\epic\\Krikwood_vs_2019_11_19_Mixed_bed_SEC.txt"), false);
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());

			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Test
	public void compareHumanToMouseEPICResult() {
		// since we want to compare components of complexes from different species, we
		// should use gene names
		ProteinComponent.TAKE_GENES_FOR_COMPARISONS = true;
		final File folder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Human_1_2_out");
		final double minOverlapScore = 0.25;
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			final File outputFile = epicComparator.compareWithOtherEPICResult(
					new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Mouse_1_2_out"),
					new File(
							"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Epic_Beta_cells_Human_VS_Mouse.txt"),
					false);
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());

			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Test
	public void compareHumanToRatEPICResult() {
		// since we want to compare components of complexes from different species, we
		// should use gene names
		ProteinComponent.TAKE_GENES_FOR_COMPARISONS = true;
		final File folder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Human_1_2_out");
		final double minOverlapScore = 0.25;
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			final File outputFile = epicComparator.compareWithOtherEPICResult(
					new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Rat_1_2_out"),
					new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Epic_Beta_cells_Human_VS_Rat.txt"),
					false);
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());

			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Test
	public void compareMouseToRatEPICResult() {
		// since we want to compare components of complexes from different species, we
		// should use gene names
		ProteinComponent.TAKE_GENES_FOR_COMPARISONS = true;
		final File folder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Mouse_1_2_out");
		final double minOverlapScore = 0.25;
		final EpicResultComparator epicComparator = new EpicResultComparator(folder, minOverlapScore);
		try {
			final File outputFile = epicComparator.compareWithOtherEPICResult(
					new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Rat_1_2_out"),
					new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Epic_Beta_cells_Mouse_VS_Rat.txt"),
					false);
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());

			System.exit(0);
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/**
	 * What is the overlapping between the proteins that form complexes in human,
	 * mouse and rat
	 */
	@Test
	public void compareHumanMouseRatEPICResults() {
		final double minOverlapScore = 0.25;
		final File humanFolder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Human_1_2_out");
		final File mouseFolder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Mouse_1_2_out");
		final File ratFolder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Beta_cell_PCP_NSAF_Rat_1_2_out");

		try {
			final File outputFile = new File(
					"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\epic\\Epic_Human_Mouse_Rat_comparison.txt");
			final FileWriter fw = new FileWriter(outputFile);
			//
			final EpicResultComparator humanEpicComparator = new EpicResultComparator(humanFolder, minOverlapScore);
			reportSummaryOfResults("Human", humanEpicComparator, fw);
			writeProteinsFile("Human", getProteins(humanEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			writeGenesFile("Human", getGenes(humanEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			//
			final EpicResultComparator mouseEpicComparator = new EpicResultComparator(mouseFolder, minOverlapScore);
			reportSummaryOfResults("Mouse", mouseEpicComparator, fw);
			writeProteinsFile("Mouse", getProteins(mouseEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			writeGenesFile("Mouse", getGenes(mouseEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			//
			final EpicResultComparator ratEpicComparator = new EpicResultComparator(ratFolder, minOverlapScore);
			reportSummaryOfResults("Rat", ratEpicComparator, fw);
			writeProteinsFile("Rat", getProteins(ratEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			writeGenesFile("Rat", getGenes(ratEpicComparator.getPredictedProteinComplexes()),
					outputFile.getParentFile());
			//
			VennData vennData = new VennData("Overlap between genes forming complexes predicted by EPIC in Beta cells",
					"Human", getGenes(humanEpicComparator.getPredictedProteinComplexes()), "Mouse",
					getGenes(mouseEpicComparator.getPredictedProteinComplexes()), "Rat",
					getGenes(ratEpicComparator.getPredictedProteinComplexes()));
			fw.write("\n" + vennData.toString() + "\n");
			//
			vennData = new VennData("Overlap between genes forming KNOWN complexes predicted by EPIC in Beta cells",
					"Human", getGenes(getKnownComplexes(humanEpicComparator.getPredictedProteinComplexes())), "Mouse",
					getGenes(getKnownComplexes(mouseEpicComparator.getPredictedProteinComplexes())), "Rat",
					getGenes(getKnownComplexes(ratEpicComparator.getPredictedProteinComplexes())));
			fw.write("\n" + vennData.toString() + "\n");
			//
			vennData = new VennData("Overlap between genes forming UNKNOWN complexes predicted by EPIC in Beta cells",
					"Human", getGenes(getUnknownComplexes(humanEpicComparator.getPredictedProteinComplexes())), "Mouse",
					getGenes(getUnknownComplexes(mouseEpicComparator.getPredictedProteinComplexes())), "Rat",
					getGenes(getUnknownComplexes(ratEpicComparator.getPredictedProteinComplexes())));
			fw.write("\n" + vennData.toString() + "\n");

			// PPIs
			final Set<String> humanPPIs = getGenesPPIs(humanEpicComparator.getProteinProteinInteractions());
			final Set<String> mousePPIs = getGenesPPIs(mouseEpicComparator.getProteinProteinInteractions());
			final Set<String> ratPPIs = getGenesPPIs(ratEpicComparator.getProteinProteinInteractions());
			vennData = new VennData(
					"Overlap between Protein-Protein Interactions forming complexes predicted by EPIC in Beta cells",
					"Human", humanPPIs, "Mouse", mousePPIs, "Rat", ratPPIs);
			fw.write("\n" + vennData.toString() + "\n");
			fw.close();
			System.out.println("Comparison report at: " + outputFile.getAbsolutePath());

			System.exit(0);
		} catch (final IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}

	private void writeGenesFile(String species, Set<String> genes, File folder) throws IOException {
		final String fileName = folder.getAbsolutePath() + File.separator + "Epic_Genes_list_" + species + ".txt";
		System.out.println("File written at: " + fileName);
		final FileWriter fw = new FileWriter(fileName);
		final List<String> list = new ArrayList<String>();
		list.addAll(genes);
		Collections.sort(list);
		for (final String gene : list) {
			fw.write(gene + "\n");
		}
		fw.close();
		System.out.println("File written at: " + fileName);
	}

	private void writeProteinsFile(String species, Set<String> proteins, File folder) throws IOException {
		final String fileName = folder.getAbsolutePath() + File.separator + "Epic_Protein_list_" + species + ".txt";
		final FileWriter fw = new FileWriter(fileName);
		final List<String> list = new ArrayList<String>();
		list.addAll(proteins);
		Collections.sort(list);
		for (final String protein : list) {
			fw.write(protein + "\n");
		}
		fw.close();
		System.out.println("File written at: " + fileName);
	}

	private Set<String> getGenesPPIs(List<ProteinProteinInteraction> proteinProteinInteractions) {
		final Set<String> ret = new THashSet<String>();
		for (final ProteinProteinInteraction ppi : proteinProteinInteractions) {
			final ProteinComponent component1 = ppi.getComponent1();
			final String gene1 = component1.getGene() != null ? component1.getGene() : component1.getKey();
			final ProteinComponent component2 = ppi.getComponent2();
			final String gene2 = component2.getGene() != null ? component2.getGene() : component2.getKey();
			final List<String> genes = new ArrayList<String>();
			genes.add(gene1);
			genes.add(gene2);
			Collections.sort(genes);
			final String genesString = genes.get(0) + "-" + genes.get(1);
			ret.add(genesString);
		}
		return ret;
	}

	private List<ProteinComplex> getKnownComplexes(List<ProteinComplex> complexes) {
		return complexes.stream().filter(c -> c.isKnown()).collect(Collectors.toList());
	}

	private List<ProteinComplex> getUnknownComplexes(List<ProteinComplex> complexes) {
		return complexes.stream().filter(c -> !c.isKnown()).collect(Collectors.toList());
	}

	private void reportSummaryOfResults(String speciesName, EpicResultComparator epicComparator, FileWriter fw)
			throws IOException {
		fw.write("\nEpic results for " + speciesName + ":\n");
		final List<ProteinComplex> usedReferenceProteinComplexes = epicComparator.getUsedReferenceProteinComplexes();
		fw.write("Number of complexes used by EPIC as reference known complexes:\t"
				+ usedReferenceProteinComplexes.size() + "\n");
		final List<ProteinProteinInteraction> ppis = epicComparator.getProteinProteinInteractions();
		fw.write("Number of protein-protein interactions:\t" + ppis.size() + "\n");
		final List<ProteinComplex> complexes = epicComparator.getPredictedProteinComplexes();
		final List<ProteinComplex> knownComplexes = getKnownComplexes(complexes);
		final List<ProteinComplex> unknownComplexes = getUnknownComplexes(complexes);
		fw.write("Number of predicted complexes:\t" + complexes.size() + "\n");
		final Set<String> proteins = getProteins(complexes);
		fw.write("Number of proteins in predicted complexes\t" + proteins.size() + "\n");
		final Set<String> genes = getGenes(complexes);
		fw.write("Number of genes in predicted complexes\t" + genes.size() + "\n");
		fw.write("Number of known predicted complexes:\t" + knownComplexes.size() + "\n");
		final Set<String> knownProteins = getProteins(knownComplexes);
		fw.write("Number of proteins in known predicted complexes\t" + knownProteins.size() + "\n");
		final Set<String> knownGenes = getGenes(knownComplexes);
		fw.write("Number of genes in known predicted complexes\t" + knownGenes.size() + "\n");
		fw.write("Number of unknown predicted complexes:\t" + unknownComplexes.size() + "\n");
		final Set<String> unknownProteins = getProteins(unknownComplexes);
		fw.write("Number of proteins in unknown predicted complexes\t" + unknownProteins.size() + "\n");
		final Set<String> unknownGenes = getGenes(unknownComplexes);
		fw.write("Number of genes in unknown predicted complexes\t" + unknownGenes.size() + "\n");

	}

	private Set<String> getProteins(List<ProteinComplex> complexes) {
		final Set<String> proteins = new THashSet<String>();
		for (final ProteinComplex complex : complexes) {
			for (final ProteinComponent component : complex.getComponentList()) {
				if (component.getAcc() != null) {
					proteins.add(component.getAcc());
				} else {
					proteins.add(component.getKey());
				}
			}
		}
		return proteins;
	}

	private Set<String> getGenes(List<ProteinComplex> complexes) {
		final Set<String> genes = new THashSet<String>();
		for (final ProteinComplex complex : complexes) {
			for (final ProteinComponent component : complex.getComponentList()) {
				if (component.getGene() != null) {
					genes.add(component.getGene());
				} else {
					genes.add(component.getKey());
				}
			}
		}
		return genes;
	}
}
