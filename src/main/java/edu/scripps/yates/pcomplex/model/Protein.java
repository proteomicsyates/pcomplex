package edu.scripps.yates.pcomplex.model;

import java.io.IOException;
import java.util.Set;

import edu.scripps.yates.pcomplex.util.DataType;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.THashSet;

public class Protein extends ProteinComponent {

	private final Set<String> others = new THashSet<String>();
	private final Double mw;
	private final TIntIntMap spcs = new TIntIntHashMap();
	private float nsaf;
	private final int fractionNumber;

	public Protein(String acc, String gene, Double mw, int spc, float nsaf, int fractionNumber) throws IOException {
		super(acc, gene);

		this.mw = mw;
		this.spcs.put(1, spc);
		this.nsaf = nsaf;
		this.fractionNumber = fractionNumber;

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
			if (p.getAcc().equals(acc) && p.getGene().equals(gene) && p.getFractionNumber() == fractionNumber) {
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
		return spcs.get(1);
	}

	public TIntIntMap getSpcs() {
		return spcs;
	}

	/**
	 * 
	 * @param rep starting by 1
	 * @return
	 */
	public Integer getSpc(int rep) {
		if (spcs.containsKey(rep)) {
			return spcs.get(rep);
		}
		return null;
	}

	public int addSpc(int spc) {
		final int repNumber = this.spcs.size() + 1;
		this.spcs.put(repNumber, spc);
		return repNumber;
	}

	public float getNSAF() {
		return nsaf;
	}

	public void setNSAF(float nsaf) {
		this.nsaf = nsaf;
	}

	public int getFractionNumber() {
		return fractionNumber;
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
