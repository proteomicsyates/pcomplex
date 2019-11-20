package edu.scripps.yates.pcomplex.model;

import java.io.IOException;
import java.util.Set;

import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.set.hash.THashSet;

public class Protein extends ProteinComponent {

	private final Set<String> others = new THashSet<String>();
	private final Double mw;
	private final int spc;
	private final float nsaf;
	private final String fractionName;

	public Protein(String acc, String gene, Double mw, int spc, float nsaf, String fractionName) throws IOException {
		super(acc, gene);

		this.mw = mw;
		this.spc = spc;
		this.nsaf = nsaf;
		this.fractionName = fractionName;

	}

	public void addOther(String otherKey) {
		others.add(otherKey);
	}

	public Set<String> getOthers() {
		return others;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Protein) {
			final Protein p = (Protein) obj;
			if (p.getAcc().equals(acc) && p.getGene().equals(gene) && p.getFractionName().equals(fractionName)) {
				if (others.isEmpty() && p.getOthers().isEmpty()) {
					return true;
				} else {
					if (others.size() != p.getOthers().size()) {
						return false;
					}
					for (final String other : others) {
						if (!p.getOthers().contains(other)) {
							return false;
						}
					}
					for (final String other : p.getOthers()) {
						if (!others.contains(other)) {
							return false;
						}
					}
					return true;
				}
			} else {
				return false;
			}
		}
		return super.equals(obj);
	}

	public Double getMw() {
		return mw;
	}

	@Override
	public String toString() {
		String string = acc + "|" + gene + "|";
		for (final String other : others) {
			string += " | " + other;
		}
		return string;
	}

	public int getSpc() {
		return spc;
	}

	public float getNSAF() {
		return nsaf;
	}

	public String getFractionName() {
		return fractionName;
	}

	public double getData(DataType dataType) {
		switch (dataType) {
		case NSAF:
			return getNSAF();
		case SPC:
			return getSpc();
		default:
			break;
		}
		return Double.NaN;
	}

}
