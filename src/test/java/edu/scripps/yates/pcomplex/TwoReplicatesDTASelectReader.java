package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.compress.utils.FileNameUtils;

import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.staticstorage.StaticProteomicsModelStorage;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Reads all dtaselects from replicate 1 and from replicate 2 and reports the
 * number of PSMs, peptides and proteins (groups) in all the data, for the 3
 * datasets (human, mouse and rat).
 * 
 * @author salvador
 *
 */
public class TwoReplicatesDTASelectReader {

//	private final static File basePath = new File(
//			"C:\\Users\\salvador\\Dropbox (Scripps Research)\\beta_cells_PCP\\DTASelect files");
	private final static File basePath = new File("D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\DTASelect files");
	private static File uniprotPath = new File("Z:\\share\\Salva\\UniprotKB");

	public static void main(String[] args) {
		try {
			final File outputFile = new File(basePath + File.separator + "Summary_numbers_1_2_replicates.txt");
			final FileWriter fw = new FileWriter(outputFile);
			final String[] species = { "Human", "Mouse", "Rat" };
			final Map<String, Set<String>> genesBySpecies = new THashMap<String, Set<String>>();
			final Map<String, Set<String>> peptidesBySpecies = new THashMap<String, Set<String>>();
			for (final String specie : species) {

				fw.write("Species:\t" + specie + "\n");
				final File[] folders = basePath.listFiles(new FileFilter() {

					@Override
					public boolean accept(File folder) {
						if (folder.isDirectory()) {
							if (FileNameUtils.getBaseName(folder.getAbsolutePath()).contains(specie)) {
								return true;
							}
						}
						return false;
					}
				});
				final List<File> dtaselectFiles = new ArrayList<File>();
				int repNum = 1;
				for (final File folder : folders) {
					final File[] dtaselectFilesInRep = folder.listFiles();
					fw.write("Fractions in replicate " + repNum + ":\t" + dtaselectFilesInRep.length + "\n");
					for (final File dtaselectFileInRep : dtaselectFilesInRep) {
						dtaselectFiles.add(dtaselectFileInRep);
					}
					repNum++;
				}

				fw.write("Total fractions:\t" + dtaselectFiles.size() + "\n");
				StaticProteomicsModelStorage.clearData();
				final DTASelectParser parser = new DTASelectParser(dtaselectFiles);
				parser.setDecoyPattern("Reverse");
				fw.write("Number of PSMs:\t" + parser.getPSMsByPSMID().size() + "\n");
				fw.write("Number of Peptides:\t" + parser.getPeptides().size() + "\n");
				final PAnalyzer panalyzer = new PAnalyzer(false);
				final List<GroupableProtein> groupableProteins = new ArrayList<GroupableProtein>();
				groupableProteins.addAll(parser.getProteins());
				final List<ProteinGroup> groups = panalyzer.run(groupableProteins);
				fw.write("Number of Protein groups:\t" + groups.size() + "\n");
				fw.write("----------------------------------------\n\n");
				fw.flush();
				// genes
				genesBySpecies.put(specie, new THashSet<String>());
				final UniprotGeneMapping geneMapper = UniprotGeneMapping.getInstance(uniprotPath, species, false, true,
						false);
				for (final Protein protein : parser.getProteins()) {
					final String acc = protein.getAccession();
					String gene = null;
					final Set<String> genes = geneMapper.mapUniprotACCToGene(acc);
					if (!genes.isEmpty()) {
						gene = genes.iterator().next();
					}
					if (gene == null) {
						gene = acc;
					}
					genesBySpecies.get(specie).add(gene);
				}
				// peptides
				peptidesBySpecies.put(specie, new THashSet<String>());
				for (final Peptide peptide : parser.getPeptides()) {
					final String sequence = peptide.getSequence();
					peptidesBySpecies.get(specie).add(sequence);
				}
			}
			reportVennDiagram("Overlap of identified genes between Human, Mouse and Rat", genesBySpecies, fw);
			reportVennDiagram("Overlap of identified peptides between Human, Mouse and Rat", peptidesBySpecies, fw);
			fw.close();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private static void reportVennDiagram(String title, Map<String, Set<String>> map, FileWriter fw)
			throws IOException {
		final List<String> species = new ArrayList<String>();
		species.addAll(map.keySet());
		final Set<String> set1 = map.get(species.get(0));
		final Set<String> set2 = map.get(species.get(1));
		final Set<String> set3 = map.get(species.get(2));
		final VennData venn = new VennData(title, species.get(0), set1, species.get(1), set2, species.get(2), set3);
		fw.write("\n\n" + venn.toString());
	}
}
