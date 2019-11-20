package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import gnu.trove.map.hash.THashMap;

public class ComplexPortalDB extends ProteinComplexDB {

	public ComplexPortalDB(File inputFile, UniprotProteinLocalRetriever uplr) throws IOException {
		super(null, "Complex Portal", false, uplr);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String line = null;
			final Map<String, Integer> indexesByHeaders = new THashMap<String, Integer>();
			int numLine = 0;
			while ((line = br.readLine()) != null) {
				numLine++;
				final String[] split = line.split("\t");
				if (numLine == 1) {
					for (int index = 0; index < split.length; index++) {
						indexesByHeaders.put(split[index], index);
					}
				} else {
					final String proteinComplexID = split[indexesByHeaders.get("#Complex ac")];
					final String proteinComplexName = split[indexesByHeaders.get("Recommended name")];
					final ProteinComplex proteinComplex = new ProteinComplex(proteinComplexID);

					proteinComplex.setName(proteinComplexName);
					final String componentsString = split[indexesByHeaders
							.get("Identifiers (and stoichiometry) of molecules in complex")];
					final String[] split2 = componentsString.split("\\|");
					for (final String component : split2) {
						// component is like Q07812(8)
						String acc = component.substring(0, component.indexOf("("));
						String gene = null;
						final Accession accession = FastaParser.getACC(acc);
						if (accession.getAccessionType() == AccessionType.UNIPROT) {
							acc = accession.getAccession();
							final Map<String, Set<String>> genes = ProteinComplexAnalyzer.getUniprotGeneMapping()
									.mapUniprotACCToGeneByType(acc);
							if (!genes.isEmpty()) {
								Set<String> genesTMP = null;
								if (genes.containsKey(UniprotGeneMapping.GENE_NAME)) {
									genesTMP = genes.get(UniprotGeneMapping.GENE_NAME);
								} else {
									genesTMP = genes.values().iterator().next();
								}
								final StringBuilder sb = new StringBuilder();
								for (final String geneTMP : genesTMP) {
									if (!"".equals(sb.toString())) {
										sb.append(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
									}
									sb.append(geneTMP);
								}
								gene = sb.toString();
							}
						} else {
							gene = acc;
							final Set<String> accs = ProteinComplexAnalyzer.getUniprotGeneMapping()
									.mapGeneToUniprotACC(gene);
							if (!accs.isEmpty()) {
								acc = accs.iterator().next();
							}

						}
						final ProteinComponent pc = new ProteinComponent(acc, gene);
						proteinComplex.addComponent(pc);
					}
					proteinComplexes.add(proteinComplex);
				}
			}
			setReady();
		} finally {
			br.close();
		}
	}

}
