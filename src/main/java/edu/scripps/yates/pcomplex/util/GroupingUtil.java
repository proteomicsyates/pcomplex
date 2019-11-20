package edu.scripps.yates.pcomplex.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import gnu.trove.set.hash.THashSet;

public class GroupingUtil {
	private final static Logger log = Logger.getLogger(GroupingUtil.class);

	public static List<ProteinGroup> grouping(File[] dtaSelectFiles, UniprotProteinLocalRetriever uplr)
			throws IOException {
		final List<File> files = new ArrayList<File>();
		for (final File file : dtaSelectFiles) {
			if (file.isFile() && file.exists()) {
				files.add(file);
			}
		}
		final DTASelectParser parser = new DTASelectParser(files);
		parser.enableProteinMergingBySecondaryAccessions(uplr, null);
		parser.setDecoyPattern("Reverse");
		final PAnalyzer panalyzer = new PAnalyzer(true);
		final List<Protein> proteins = parser.getProteins();
		log.info(proteins.size() + " proteins read from " + files.size() + " files. Now grouping them with PANalyzer");
		final Set<GroupableProtein> groupableProteins = new THashSet<GroupableProtein>();
		for (final GroupableProtein groupableProtein : proteins) {
			groupableProteins.add(groupableProtein);
		}
		final List<ProteinGroup> proteinGroups = panalyzer.run(groupableProteins);
		log.info(proteinGroups.size() + " protein groups found");
		final List<ProteinGroup> ret = proteinGroups.stream()
				.filter(group -> group.getEvidence() != ProteinEvidence.NONCONCLUSIVE).collect(Collectors.toList());
		log.info(proteinGroups.size() + " protein groups after removing non conclusive proteins");
		return ret;

	}

}
