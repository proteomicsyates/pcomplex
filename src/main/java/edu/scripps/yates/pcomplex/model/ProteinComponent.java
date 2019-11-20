package edu.scripps.yates.pcomplex.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.builder.HashCodeBuilder;

import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;

public class ProteinComponent {
	protected String acc;
	protected String gene;
	private String proteinName;
	private String toString;
	private int hashCode = -1;
	// to use such as when there is an ambiguity, we remove one by keeping the
	// protein that is in the reference database.
	public static List<ProteinComplexDB> dbs = new ArrayList<ProteinComplexDB>();

	public ProteinComponent(String acc, String gene) throws IOException {
		this.acc = chooseOne(acc);
		this.gene = chooseOne(gene);
		// deal with ambiguities and select one
		if (acc == null || gene == null) {
			if (gene != null) {
				final Set<String> accs = ProteinComplexAnalyzer.getUniprotGeneMapping().mapGeneToUniprotACC(gene);
				if (!accs.isEmpty()) {
					for (final String acc2 : accs) {
						if (isInReferenceDB(acc)) {
							this.acc = acc2;
						}
					}
					if (this.acc == null) {
						this.acc = accs.iterator().next();
					}
				} else {
					this.acc = null;
				}
				this.gene = gene;
			} else {
				final Set<String> genes = ProteinComplexAnalyzer.getUniprotGeneMapping().mapUniprotACCToGene(acc);
				if (!genes.isEmpty()) {
					for (final String gene2 : genes) {
						if (isInReferenceDB(acc)) {
							this.gene = gene2;
						}
					}
					if (this.gene == null) {
						this.gene = genes.iterator().next();
					}
				} else {
					this.gene = null;
				}
				this.acc = acc;
			}
		} else {
			this.acc = acc;
			this.gene = gene;
		}

	}

	private boolean isInReferenceDB(String acc) {

		for (final ProteinComplexDB db : dbs) {
			if (!db.getProteinComplexesByProtein(acc).isEmpty()) {
				return true;
			}
		}

		return false;
	}

	/**
	 * If there is more than one acc, Returns the one that is present in the
	 * reference database or the first one if not
	 * 
	 * @param acc2
	 * @return
	 */
	private String chooseOne(String acc2) {
		if (acc2 != null && acc2.contains(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR)) {
			final String[] split = acc2.split(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
			for (final String acc : split) {
				if (isInReferenceDB(acc)) {
					return acc;
				}
			}
			return split[0];
		} else {
			return acc2;
		}
	}

	public String getAcc() {
		return acc;
	}

	public String getGene() {
		return gene;
	}

	@Override
	public int hashCode() {
		if (hashCode == -1) {
			hashCode = HashCodeBuilder.reflectionHashCode(getKey());
		}
		return hashCode;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ProteinComponent) {
			// 2 components are true if they have the same ACC or Gene, taking
			// into account that they can be ambiguous
			final ProteinComponent pc = (ProteinComponent) obj;
			return pc.getKey().equals(getKey());
		}
		return super.equals(obj);
	}

	@Override
	public String toString() {
		if (toString == null) {
			toString = acc + "|" + gene;
		}
		return toString;
	}

	public String toString(boolean useGeneToPrint) {
		if (useGeneToPrint) {
			return gene;
		} else {
			return acc;
		}
	}

	public String getKey() {
		if (acc != null) {
			return acc;
		}
		return gene;
	}

	public String getProteinName() {
		return proteinName;
	}

	public void setProteinName(String proteinName) {
		this.proteinName = proteinName;
	}

}
