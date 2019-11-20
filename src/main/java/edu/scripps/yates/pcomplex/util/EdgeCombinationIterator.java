package edu.scripps.yates.pcomplex.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.CombinatoricsUtils;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.set.hash.THashSet;
import smile.math.Math;

public class EdgeCombinationIterator implements Iterator<List<Edge>> {
	private final List<Edge> edges = new ArrayList<Edge>();
	private final Iterator<int[]> combinationsIterator;
	private final ProgressCounter progress;

	public EdgeCombinationIterator(int size, Collection<Edge> edges) {
		this.edges.addAll(edges);
		combinationsIterator = CombinatoricsUtils.combinationsIterator(edges.size(), size);
		Double total = Math.factorial(edges.size()) / (Math.factorial(size) * Math.factorial(edges.size() - size));
		if (Double.isNaN(total)) {
			total = getNumCombinations(edges.size(), size);
		}
		if (!Double.isNaN(total)) {
			progress = new ProgressCounter(total.longValue(), ProgressPrintingType.PERCENTAGE_STEPS, 0);
			progress.setShowRemainingTime(true);
			progress.setSuffix("Combinations generated for MMR calculation");
		} else {
			progress = null;
		}

	}

	public static void main(String[] args) {
		final double num = getNumCombinations(752, 14);
		System.out.println(num);
		System.out.println(Double.valueOf(num).longValue());
	}

	public static double getNumCombinations(int n, int k) {
		double ret = 1.0;

		for (int i = 0; i < k; i++) {
			ret = ret * ((n - i) / (i + 1));
			if (ret < 0.0) {
				throw new IllegalArgumentException("double overflow");
			}
		}
		return ret;
	}

	@Override
	public boolean hasNext() {

		return combinationsIterator.hasNext();
	}

	public String getProgress() {
		if (progress == null) {
			return "";
		}
		return progress.printIfNecessary();
	}

	@Override
	public List<Edge> next() {
		final int[] combinationsIndexes = combinationsIterator.next();
		if (progress != null) {
			progress.increment();
		}
		final List<Edge> list = new ArrayList<Edge>();
		final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
		for (final int index : combinationsIndexes) {
			final Edge edge = edges.get(index);
			if (set.contains(edge.getPredicted()) || set.contains(edge.getReference())) {
				return null;
			}
			set.add(edge.getPredicted());
			set.add(edge.getReference());
			list.add(edge);
		}
		return list;

	}

}
