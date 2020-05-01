package edu.scripps.yates.pcomplex.model;

public class ProteinProteinInteraction {
	private final ProteinComponent component1;
	private final ProteinComponent component2;
	private final double score;

	public ProteinProteinInteraction(ProteinComponent component1, ProteinComponent component2, double score) {
		super();
		this.component1 = component1;
		this.component2 = component2;
		this.score = score;
	}

	public ProteinComponent getComponent1() {
		return component1;
	}

	public ProteinComponent getComponent2() {
		return component2;
	}

	public double getScore() {
		return score;
	}

}
