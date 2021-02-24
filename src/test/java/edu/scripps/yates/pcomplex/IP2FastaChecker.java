package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import com.compomics.dbtoolkit.io.implementations.FASTADBLoader;
import com.compomics.util.protein.Protein;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.SftpException;

import edu.scripps.yates.pcomplex.ftp.MySftpProgressMonitor;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.set.hash.THashSet;

public class IP2FastaChecker {
	private final File propertiesFile = new File("Z:\\share\\Salva\\data\\IP2_Fastas\\ip2FileTransfer.properties");
	private final IP2Util ip2Util = new IP2Util(new MySftpProgressMonitor(System.out), propertiesFile, "");
	private final File databaseListCVSFile = new File("Z:\\share\\Salva\\data\\IP2_Fastas\\Database_list.csv");

	@Test
	public void test() {
		try {
			final FileWriter fw = new FileWriter("Z:\\share\\Salva\\data\\IP2_Fastas\\report.txt");
			final Set<String> users = getUsers();
			for (final String user : users) {
				final List<String> fastaFiles = new ArrayList<String>();

				fastaFiles.addAll(ip2Util.getFastaFiles(user));
				System.out.println("User " + user + " has " + fastaFiles.size() + " databases");
				fw.write("user\tFASTA\ttotal\tnum forward\tnum decoy\tdifference\terror_reading\n");
				for (final String fastaFile : fastaFiles) {
					final StringBuilder report = new StringBuilder();

					final File outputFile = new File(
							"Z:\\share\\Salva\\data\\IP2_Fastas" + File.separator + FilenameUtils.getName(fastaFile));
					if (!outputFile.exists()) {
						outputFile.getParentFile().mkdirs();
						final OutputStream outputStream = new FileOutputStream(outputFile, false);
						ip2Util.download(fastaFile, outputStream);
						System.out.println(outputFile.getAbsolutePath() + " created");
					}
					System.out.println("Reading fasta " + FilenameUtils.getName(fastaFile));
					final FASTADBLoader loader = new FASTADBLoader();
					try {
						report.append(user + "\t" + FilenameUtils.getName(fastaFile) + "\t");
						loader.load(outputFile.getAbsolutePath());
						Protein protein;
						int numFW = 0;
						int numREV = 0;
						while ((protein = loader.nextProtein()) != null) {
							final boolean reverse = FastaParser.isReverse(protein.getHeader().getRawHeader());
							if (reverse) {
								numREV++;
							} else {
								numFW++;
							}
						}
						final int total = numFW + numREV;

						report.append(total + "\t" + numFW + "\t" + numREV);
						if (numFW != numREV) {
							report.append("\t" + (numFW - numREV));
						}
					} catch (final Exception e) {
						e.printStackTrace();
						report.append("\t\t\t\t" + e.getMessage());
					}
					report.append("\n");
					fw.write(report.toString());
					fw.flush();
					System.out.println(report.toString());
				}
			}
			fw.close();
		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final JSchException e) {
			e.printStackTrace();
		} catch (final SftpException e) {
			e.printStackTrace();
		}
	}

	private Set<String> getUsers() throws IOException {
		final List<String> users = FileUtils.readColumnFromTextFile(databaseListCVSFile, ",", 0, false);
		final Set<String> userSet = new THashSet<String>();
		userSet.addAll(users);
		System.out.println(userSet.size() + " users");
		return userSet;
	}
}
