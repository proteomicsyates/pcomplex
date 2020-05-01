package edu.scripps.yates.pcomplex.gprofiler;

import java.util.Set;

import gnu.trove.set.hash.THashSet;

public class GOTerm {
	private final Set<String> accs = new THashSet<String>();
	private final GO_TYPE goType;
	private final String termName;
	private final String goID;
	private final double adjPValue;
	private final int termSize;
	private final int querySize;
	private final int intersectionSize;
	private final int effectiveDomainSize;

	public GOTerm(String[] gProfilerSplittedLine) {
		this(gProfilerSplittedLine[0], gProfilerSplittedLine[1], gProfilerSplittedLine[2], gProfilerSplittedLine[3],
				gProfilerSplittedLine[5], gProfilerSplittedLine[6], gProfilerSplittedLine[7], gProfilerSplittedLine[8],
				gProfilerSplittedLine[9]);
	}

	public GOTerm(String goTypeString, String termName, String goID, String adjPValue, String termSize,
			String querySize, String intersectionSize, String effectiveDomainSize, String intersectionsString) {
		super();
		if (goTypeString.startsWith("\"")) {
			goTypeString = goTypeString.substring(1);
		}
		this.goType = GO_TYPE.parseString(goTypeString);
		this.termName = termName;
		this.goID = goID;
		this.adjPValue = Double.valueOf(adjPValue);
		this.termSize = Integer.valueOf(termSize);
		this.querySize = Integer.valueOf(querySize);
		this.intersectionSize = Integer.valueOf(intersectionSize);
		this.effectiveDomainSize = Integer.valueOf(effectiveDomainSize);
		if (intersectionsString.contains(",")) {
			final String[] split = intersectionsString.split(",");
			for (String acc : split) {
				if (acc.endsWith("\"")) {
					acc = acc.substring(0, acc.length() - 1);
				}
				addAcc(acc);
			}
		} else {
			if (intersectionsString.endsWith("\"")) {
				intersectionsString = intersectionsString.substring(0, intersectionsString.length() - 1);
			}
			addAcc(intersectionsString);
		}
	}

	public boolean containsProtein(String acc) {
		return getAccs().contains(acc);
	}

	private void addAcc(String acc) {
		this.accs.add(acc);
	}

	public Set<String> getAccs() {
		return accs;
	}

	public GO_TYPE getGoType() {
		return goType;
	}

	public String getTermName() {
		return termName;
	}

	public String getGoID() {
		return goID;
	}

	public double getAdjPValue() {
		return adjPValue;
	}

	public int getTermSize() {
		return termSize;
	}

	public int getQuerySize() {
		return querySize;
	}

	public int getIntersectionSize() {
		return intersectionSize;
	}

	public int getEffectiveDomainSize() {
		return effectiveDomainSize;
	}

}
