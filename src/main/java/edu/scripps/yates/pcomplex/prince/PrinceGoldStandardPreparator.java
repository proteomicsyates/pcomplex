package edu.scripps.yates.pcomplex.prince;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

import org.apache.commons.io.FilenameUtils;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.db.ComplexPortalDB;

/**
 * Prince needs a gold standard for its analysis which is know protein
 * complexes. This tool will read the Core, ComplexPortal and QuickGO to
 * elaborate a file with the complexes. This file will be read by a script in R
 * to be incorporated in Prince analysis in R.
 * 
 * @author salvador
 *
 */
public class PrinceGoldStandardPreparator {

	public static void main(String[] args) {
		final PrinceGoldStandardPreparator p = new PrinceGoldStandardPreparator();
		try {
			p.run();
			System.out.println("Everything ok");
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private final File uniprotFolder = new File("Z:\\share\\Salva\\UniprotKB");

	private final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotFolder, true);

	private void run() throws IOException {
		final File complexPortalFolder = new File(
				"D:\\Dropbox (Scripps Research)\\beta_cells_PCP\\databases\\EBI-Complex-Portal");
		processEBIComplexPortal(complexPortalFolder);

	}

	private void processEBIComplexPortal(File complexPortalFolder) throws IOException {
		final File[] complexFiles = getComplexPortalFiles(complexPortalFolder);
		for (final File complexFile : complexFiles) {
			final ComplexPortalDB db = new ComplexPortalDB(complexFile, uplr);
			final File outputFile = new File(complexFile.getParent() + File.separator
					+ FilenameUtils.getBaseName(complexFile.getAbsolutePath()) + "_prince.tsv");
			db.exportComplexes(outputFile);

		}

	}

	private File[] getComplexPortalFiles(File complexPortalFolder) {
		return complexPortalFolder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.contains("prince")) {
					return false;
				}
				return true;
			}
		});

	}
}
