package edu.scripps.yates.pcomplex;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.ChannelSftp.LsEntry;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpException;

import edu.scripps.yates.pcomplex.ftp.MySftpProgressMonitor;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.ftp.FTPUtils;
import edu.scripps.yates.utilities.properties.PropertiesUtil;
import gnu.trove.map.hash.THashMap;

public class IP2Util {
	private final static Logger log = Logger.getLogger(IP2Util.class);
	private final MySftpProgressMonitor progressMonitor;
	private final File propertiesFile;
	private final String projectName;
	private final String projectFullPath;
	public final static String IP2_SERVER_PROJECT_BASE_PATH = "ip2_server_project_base_path";
	public final static String PROJECT_NAME = "project_name";

	public IP2Util(MySftpProgressMonitor progressMonitor, File propertiesFile, String projectName) {

		this.progressMonitor = progressMonitor;
		this.propertiesFile = propertiesFile;
		this.projectName = projectName.replace(" ", "_");
		projectFullPath = getProperties(this.propertiesFile).getProperty(IP2_SERVER_PROJECT_BASE_PATH) + "/"
				+ this.projectName;

		// try {
		// final List<Project> allProjects =
		// ProjectService.getAllProjects("ip2");
		// for (final Project project : allProjects) {
		// final Integer projectID = project.getId();
		// final List<MspExperiment> allExperiments =
		// ExperimentService.getAllExperiments(projectID);
		// for (final MspExperiment experiment : allExperiments) {
		// if (experiment.getId() == 242600) {
		// final List<DbSearch> dbSearches = experiment.getDbSearch();
		// for (final DbSearch dbSearch : dbSearches) {
		// log.info(dbSearch);
		// }
		// }
		// }
		// }
		// } catch (final APIException e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// }
	}

	protected Session loginToIP2() throws IOException {
		log.info("Login into server...");
		final Properties properties = getProperties(propertiesFile);
		final String hostName = properties.getProperty("ip2_server_url");
		final String userName = properties.getProperty("ip2_server_user_name");
		final String password = properties.getProperty("ip2_server_password");
		final int port = Integer.valueOf(properties.getProperty("ip2_server_connection_port"));
		try {
			final Session ftpClient = FTPUtils.loginSSHClient(hostName, userName, password, port);
			log.info("Login OK");
			return ftpClient;
		} catch (final JSchException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Cannot connect to IP2: " + e.getMessage());
		}

	}

	protected static Properties getProperties(File propertiesFile) {
		try {
			final Properties properties = PropertiesUtil.getProperties(propertiesFile);
			return properties;
		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e);
			throw new IllegalArgumentException("Error reading properties file: " + propertiesFile.getAbsolutePath());
		}
	}

	public Map<String, String> getDTASelectsInProject() throws JSchException, IOException, SftpException {
		final Map<String, String> ret = new THashMap<String, String>();
		Session sftpIP2 = null;
		try {
			sftpIP2 = loginToIP2();

			final ChannelSftp sftpChannel = FTPUtils.openSFTPChannel(sftpIP2);
			final Vector<LsEntry> ls = sftpChannel.ls(projectFullPath);
			final List<String> experimentNames = new ArrayList<String>();
			for (final LsEntry lsEntry : ls) {
				final String experimentName = lsEntry.getFilename();
				if (experimentName.contains("HEK") || experimentName.contains("hplc") || experimentName.contains("AS")
						|| experimentName.contains("AJS") || experimentName.startsWith("Cond_")
						|| experimentName.startsWith("PT") || experimentName.startsWith("AJS_MB_SEC")
						|| experimentName.startsWith("AJS_MB")) {
					experimentNames.add(experimentName);
				} else {
					continue;
				}
			}
			if (experimentNames.isEmpty()) {
				throw new IllegalArgumentException(
						"Experiment names are not recognized. Come just above this line and modify it accordingly");
			}
			for (final String experimentName : experimentNames) {

				// final int frac =
				// getFractionNumberFromExperimentName(experimentName);
				final String searchesPath = projectFullPath + "/" + experimentName + "/search";
				final Vector<LsEntry> ls2 = sftpChannel.ls(searchesPath);
				for (final LsEntry lsEntry : ls2) {
					final String name = lsEntry.getFilename();
					if (!name.equals(".") && !name.equals("..")) {
						final String searchPath = projectFullPath + "/" + experimentName + "/search/" + name;
						final Vector<LsEntry> ls3 = sftpChannel.ls(searchPath);
						for (final LsEntry lsEntry2 : ls3) {
							final String name2 = lsEntry2.getFilename();
							if (name2.equals("DTASelect-filter.txt")) {
								// if the experiment already have a search, take
								// the newest search
								if (ret.containsKey(experimentName)) {
									String searchPath1 = ret.get(experimentName);
									searchPath1 = searchPath1.substring(0,
											searchPath1.indexOf("/DTASelect-filter.txt"));
									final int searchNumber1 = Integer
											.valueOf(searchPath1.substring(searchPath1.lastIndexOf("_") + 1));
									final int searchNumber2 = Integer
											.valueOf(searchPath.substring(searchPath.lastIndexOf("_") + 1));
									if (searchNumber2 > searchNumber1) {
										ret.put(experimentName, searchPath + "/DTASelect-filter.txt");
										log.info("Avoiding to override newer search in " + searchesPath + " where "
												+ searchNumber2 + " is newer than " + searchNumber1);
									} else {
										log.info("Avoiding to override newer search in " + searchesPath + " where "
												+ searchNumber1 + " is newer than " + searchNumber2);
									}
								} else {
									ret.put(experimentName, searchPath + "/DTASelect-filter.txt");
								}
							}
						}
					}
				}
			}

			// sftpChannel.get(fullPathInIP2, outputStream);
			sftpChannel.exit();

			log.info(ret.size() + " different dtaselects out of " + experimentNames.size() + "experiment names");

		} finally {
			if (sftpIP2 != null) {
				sftpIP2.disconnect();
			}
		}
		return ret;
	}

	public static int getFractionNumberFromExperimentName(String experimentName, String experimentNamePattern) {
		try {
			if (experimentName.startsWith("HEK293_WCX_")) {
				final String tmp = experimentName.substring(11);
				String tmp2 = tmp.substring(0, tmp.indexOf("_"));
				if (tmp2.contains("frac")) {
					tmp2 = tmp2.replace("frac", "");
				}
				return Integer.valueOf(tmp2);
			} else if (experimentName.startsWith("HEK293_")) {
				final String tmp = experimentName.substring(7);
				String tmp2 = tmp.substring(0, tmp.indexOf("_"));
				if (tmp2.contains("frac")) {
					tmp2 = tmp2.replace("frac", "");
				}
				return Integer.valueOf(tmp2);
			} else if (experimentName.contains("AS_fraction_")) {
				final String tmp = experimentName.substring(experimentName.indexOf("AS_fraction_") + 12);
				final String tmp2 = tmp.substring(0, tmp.indexOf("_"));
				return Integer.valueOf(tmp2);
			} else if (experimentName.contains("frac")) {

				final String tmp = experimentName.substring(experimentName.indexOf("frac") + 4);
				final String tmp2 = tmp.substring(0, tmp.indexOf("_"));
				return Integer.valueOf(tmp2);
			} else if (experimentName.contains(experimentNamePattern)) {

				final String tmp = experimentName.substring(experimentName.indexOf("Frac") + 5);
				final String tmp2 = tmp.substring(0, tmp.indexOf("_"));
				return Integer.valueOf(tmp2);
			} else {
				// such as 20190123_AJS_SEC-24
				final Pattern patter = Pattern.compile(".*SEC_(\\d+)_.*");
				final Matcher matcher = patter.matcher(experimentName);
				if (matcher.find()) {
					final String tmp = matcher.group(1);
					return Integer.valueOf(tmp);
				}
				// such as 20190123_AJS_WCX_45_2019_02_01_14_242676
				final Pattern patter2 = Pattern.compile(".*WCX_(\\d+)_.*");
				final Matcher matcher2 = patter2.matcher(experimentName);
				if (matcher2.find()) {
					final String tmp = matcher2.group(1);
					return Integer.valueOf(tmp);
				}
				// such as AJS_MB_SEC_100
				final Pattern patter3 = Pattern.compile(".*SEC_(\\d+)");
				final Matcher matcher3 = patter3.matcher(experimentName);
				if (matcher3.find()) {
					final String tmp = matcher3.group(1);
					return Integer.valueOf(tmp);
				}
				// such as PT1781S1F16
				final Pattern patter4 = Pattern.compile("\\w+F(\\d+)");
				final Matcher matcher4 = patter4.matcher(experimentName);
				if (matcher4.find()) {
					final String tmp = matcher4.group(1);
					return Integer.valueOf(tmp);
				}
				// such as AJS_MB_28_2020_02_25_10_252010
				final Pattern patter5 = Pattern.compile("AJS_MB_(\\d+)_.*");
				final Matcher matcher5 = patter5.matcher(experimentName);
				if (matcher5.find()) {
					final String tmp = matcher5.group(1);
					return Integer.valueOf(tmp);
				}
				throw new IllegalArgumentException(
						"It was not possible to extract fraction number from experiment name: '" + experimentName
								+ "'");
			}
		} catch (final NumberFormatException e) {
			throw new IllegalArgumentException(
					"It was not possible to extract fraction number from experiment name: '" + experimentName + "'");
		}
	}

	public void download(String fullPathToIP2, OutputStream outputStream)
			throws IOException, JSchException, SftpException {
		final Session sftpIP2 = loginToIP2();
		final long bytes = FTPUtils.download(sftpIP2, fullPathToIP2, outputStream, progressMonitor);
		log.info("File downloaded " + FileUtils.getDescriptiveSizeFromBytes(bytes));
		sftpIP2.disconnect();
	}

	public List<String> getFastaFiles(String user) throws IOException, JSchException, SftpException {
		final List<String> ret = new ArrayList<String>();
		Session sftpIP2 = null;
		try {
			sftpIP2 = loginToIP2();

			final ChannelSftp sftpChannel = FTPUtils.openSFTPChannel(sftpIP2);

			Vector<LsEntry> ls = null;
			final String path = projectFullPath + user + "/database";
			try {

				ls = sftpChannel.ls(path);
				for (final LsEntry lsEntry : ls) {
					final String fastaFile = lsEntry.getFilename();
					if (fastaFile.endsWith(".fasta")) {
						ret.add(path + "/" + fastaFile);
					} else {
						continue;
					}
				}
			} catch (final Exception e) {
				log.warn(e.getMessage());
				log.warn(user + " has no databases at " + path);
			}
			// sftpChannel.get(fullPathInIP2, outputStream);
			sftpChannel.exit();
			log.info("Transfer done. " + ret.size() + " fastas");

		} finally {
			if (sftpIP2 != null) {
				sftpIP2.disconnect();
			}
		}
		return ret;
	}
}
