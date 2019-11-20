package edu.scripps.yates.pcomplex.util;

import edu.scripps.yates.pcomplex.model.ProteinComplex;

public class Edge {
	private final ProteinComplex reference;
	private final ProteinComplex predicted;
	private final double weight;

	Edge(ProteinComplex reference, ProteinComplex predicted, double weight) {
		this.predicted = predicted;
		this.reference = reference;
		this.weight = weight;
	}

	public double getWeight() {
		return weight;
	}

	public ProteinComplex getReference() {
		return reference;
	}

	public ProteinComplex getPredicted() {
		return predicted;
	}

	@Override
	public String toString() {
		return "Edge[" + reference + "," + predicted + ", w=" + weight + "]";
	}
}
