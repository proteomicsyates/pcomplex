package edu.scripps.yates.pcomplex.dbs;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.db.IntActDB;

public class IntActTest {

	@Test
	public void intActTest() {

		final File uniprotFolder = new File("C:\\Users\\salvador\\Desktop\\uniprotKB");
		try {
			final String[] taxIDs = { "9606", "10090", "10116" };
			final String[] speciesList = { "homo_sapiens", "mus_musculus", "rattus_norvegicus" };
			for (int i = 0; i < taxIDs.length; i++) {
				final String taxID = taxIDs[i];
				final String species = speciesList[i];
				final File intActFile = new File(
						"C:\\Users\\salvador\\Dropbox (Scripps Research)\\beta_cells_PCP\\databases\\IntAct\\IntAct_"
								+ species + ".tsv");
				final IntActDB dbHuman = new IntActDB(intActFile, taxID,
						new UniprotProteinLocalRetriever(uniprotFolder, true));
				System.out.println("*********************************");
				System.out.println(dbHuman.getProteinComplexes().size() + " complexes in " + taxID);
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}
}
