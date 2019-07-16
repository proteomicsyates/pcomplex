package edu.scripps.yates.pcomplex.cofractionation;

/**
 * Types of complexes, INTRA_COMPLEX, INTER_COMPLEX or NOVEL INTERACTORS
 * 
 * @author salvador
 *
 */
public enum ClassLabel {
	INTRA_COMPLEX(0, "protein pairs that occur in the same complex"), INTER_COMPLEX(1,
			"Proteins that are known to interact but not between them"), NOVEL_INTERACTORS(2,
					"One or both proteins are not known to interact");
	private final int classNumber;
	private final String explanation;

	private ClassLabel(int classNumber, String explanation) {
		this.classNumber = classNumber;
		this.explanation = explanation;
	}

	public int getClassNumber() {
		return classNumber;
	}

	public String getExplanation() {
		return explanation;
	}

	@Override
	public String toString() {
		return getExplanation();
	}
}
