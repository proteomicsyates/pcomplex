package edu.scripps.yates.pcomplex.gprofiler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class GProfilerResult {
	private final List<GOTerm> goTerms = new ArrayList<GOTerm>();
	private final Map<String, List<GOTerm>> goTermsByProteinAcc = new THashMap<String, List<GOTerm>>();
	private final Map<String, Set<String>> proteinsByGOID = new THashMap<String, Set<String>>();
	private final Map<GO_TYPE, List<GOTerm>> goTermsByGOType = new THashMap<GO_TYPE, List<GOTerm>>();

	public void addGOTerm(GOTerm go) {
		goTerms.add(go);
		final String goID = go.getGoID();
		if (!proteinsByGOID.containsKey(goID)) {
			proteinsByGOID.put(goID, new THashSet<String>());
		}
		for (final String acc : go.getAccs()) {
			if (!goTermsByProteinAcc.containsKey(acc)) {
				goTermsByProteinAcc.put(acc, new ArrayList<GOTerm>());
			}
			goTermsByProteinAcc.get(acc).add(go);
			proteinsByGOID.get(goID).add(acc);
		}
		if (!goTermsByGOType.containsKey(go.getGoType())) {
			goTermsByGOType.put(go.getGoType(), new ArrayList<GOTerm>());
		}
		goTermsByGOType.get(go.getGoType()).add(go);

	}

	public List<GOTerm> getGOTermsByProtein(String acc) {
		final List<GOTerm> ret = goTermsByProteinAcc.get(acc);
		if (ret != null) {
			return ret;
		} else {
			return Collections.emptyList();
		}
	}

	public int getNumGOTerms() {
		return goTerms.size();
	}

	public List<GOTerm> getGOTermsByGOType(GO_TYPE goType) {
		return goTermsByGOType.get(goType);
	}

	public Set<String> getProteinsWithGOTerm(String goID) {
		final Set<String> set = proteinsByGOID.get(goID);
		if (set != null) {
			return set;
		}
		return Collections.emptySet();

	}
}
