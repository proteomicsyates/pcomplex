package edu.scripps.yates.pcomplex.epic;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.model.ProteinProteinInteraction;
import gnu.trove.set.hash.THashSet;

public class EpicPPIReader {
	private final static Logger log = Logger.getLogger(EpicPPIReader.class);

	public static List<ProteinProteinInteraction> readPPIs(File ppisEPICFile) throws IOException {
		return EpicResultComparator.readProteinProteinInteractionsFromExpClustFile(ppisEPICFile);
	}

	public static void main(String[] args) {
		try {
			final List<ProteinProteinInteraction> ppis = EpicPPIReader
					.readPPIs(new File(args[0] + File.separator + "Out.rf.exp.pred.txt"));
			final File output = new File(
					args[0] + File.separator + FilenameUtils.getBaseName(args[0]) + "_ppi_separated.txt");
			final Set<String> accs = new THashSet<String>();
			for (final ProteinProteinInteraction ppi : ppis) {
				accs.addAll(ppi.getComponent1().getInvidualAccs());
				accs.addAll(ppi.getComponent2().getInvidualAccs());
			}
			final FileWriter fw = new FileWriter(output);
			for (final String acc : accs) {
				fw.write(acc + "\n");
			}
			fw.close();
			log.info(accs.size() + " proteins from the interactions written in file " + output.getAbsolutePath());

			// GO analyzer

		} catch (final IOException e) {
			e.printStackTrace();
		}
	}
}
