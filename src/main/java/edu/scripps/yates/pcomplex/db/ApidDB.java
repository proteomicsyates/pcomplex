package edu.scripps.yates.pcomplex.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class ApidDB extends ProteinComplexDB {

	private final static Logger log = Logger.getLogger(ApidDB.class);
	private final String complexesFileName = "_complexes.txt";
	private final File inputFile;

	public ApidDB(File inputFile, UniprotProteinLocalRetriever uplr) throws IOException {
		super(null, "APIDdb", false, uplr, true);
		this.inputFile = inputFile;
		final File complexesFile = new File(inputFile.getParentFile().getAbsolutePath() + File.separator
				+ FilenameUtils.getBaseName(inputFile.getAbsolutePath()) + complexesFileName);
		if (!complexesFile.exists() || complexesFile.length() == 0l) {
			BufferedReader br = null;
			try {

				br = new BufferedReader(new FileReader(inputFile));
				String line = null;
				int numLine = 1;
				TObjectIntMap<String> indexesByHeader = new TObjectIntHashMap<String>();
				while ((line = br.readLine()) != null) {
					if (numLine == 1) {
						indexesByHeader = getIndexesByHeader(line);
						numLine++;
					} else {
						if (line.contains("\t")) {
							final String[] split = line.split("\t");
							final ProteinComplex proteinComplex = new ProteinComplex(String.valueOf(numLine++));
							final String uniprotA = split[indexesByHeader.get("UniprotID_A")];
							final String geneA = split[indexesByHeader.get("GeneName_A")];
							final ProteinComponent componentA = new ProteinComponent(uniprotA, geneA);

							final String uniprotB = split[indexesByHeader.get("UniprotID_B")];
							final String geneB = split[indexesByHeader.get("GeneName_B")];
							final ProteinComponent componentB = new ProteinComponent(uniprotB, geneB);
							if (componentA.equals(componentB)) {
								continue;
							}
							proteinComplex.addComponent(componentA);
							proteinComplex.addComponent(componentB);
							proteinComplexes.add(proteinComplex);
						}
					}
				}
				setReady();
			} finally {
				br.close();
			}
		} else {
			BufferedReader br = null;
			try {

				br = new BufferedReader(new FileReader(complexesFile));
				String line = null;
				while ((line = br.readLine()) != null) {

					final String[] split = line.split("\t");
					final ProteinComplex proteinComplex = new ProteinComplex(split[0]);
					for (int i = 1; i < split.length; i++) {
						final String accGene = split[i];
						final String[] split2 = accGene.split("|");
						final String acc = split2[0];
						final String gene = split2[1];

						final ProteinComponent pc = new ProteinComponent(acc, gene);
						proteinComplex.addComponent(pc);

					}
					proteinComplexes.add(proteinComplex);
				}
				setIndexByGene(true);
				super.setReady();
			} finally {
				br.close();
			}
		}
	}

	private TObjectIntMap<String> getIndexesByHeader(String line) {
		final TObjectIntMap<String> ret = new TObjectIntHashMap<String>();
		final String[] split = line.split("\t");
		for (int index = 0; index < split.length; index++) {
			ret.put(split[index], index);
		}
		return ret;
	}

	@Override
	public void setReady() {

		super.setReady();
		// write to complexes file
		writeToComplexesFile();
	}

	private void writeToComplexesFile() {
		final File complexesFile = new File(inputFile.getParentFile().getAbsolutePath() + File.separator
				+ FilenameUtils.getBaseName(inputFile.getAbsolutePath()) + complexesFileName);
		FileWriter fw = null;
		try {
			fw = new FileWriter(complexesFile);
			for (final ProteinComplex proteinComplex : proteinComplexes) {
				fw.write(proteinComplex.getId() + "\t");
				final StringBuilder sb = new StringBuilder();
				for (final ProteinComponent component : proteinComplex.getComponentList()) {
					if (!"".equals(sb.toString())) {
						sb.append("\t");
					}
					sb.append(component.getAcc() + "|" + component.getGene());
				}
				fw.write(sb.toString());
				fw.write("\n");
			}

		} catch (final IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error writting complex file: " + complexesFile.getAbsolutePath());
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (final IOException e) {
					e.printStackTrace();
				}
			}
		}

	}

}
