package edu.scripps.yates.pcomplex.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ClusterEvaluation {
	private final static Logger log = Logger.getLogger(ClusterEvaluation.class);

	public static double getOverlap(ProteinComplex complex1, ProteinComplex complex2) {
		final double overlap = getOverlap(complex1.getComponentList(), complex2.getComponentList());
		return overlap;
	}

	public static Set<ProteinComponent> getIntersection(ProteinComplex complex1, ProteinComplex complex2) {
		return getIntersection(complex1.getComponentList(), complex2.getComponentList());
	}

	public static Set<ProteinComponent> getIntersection(List<ProteinComponent> componentList1,
			List<ProteinComponent> componentList2) {
		final Set<ProteinComponent> intersection = new THashSet<ProteinComponent>();

		List<ProteinComponent> smallComplex = componentList1;
		List<ProteinComponent> bigComplex = componentList2;
		if (smallComplex.size() > bigComplex.size()) {
			smallComplex = componentList2;
			bigComplex = componentList1;
		}
		final Set<ProteinComponent> smallSet = new THashSet<ProteinComponent>();
		smallSet.addAll(smallComplex);
		for (final ProteinComponent c1 : bigComplex) {
			if (smallSet.contains(c1)) {
				intersection.add(c1);
			}
		}

		return intersection;
	}

	private static double getOverlap(List<ProteinComponent> componentList1, List<ProteinComponent> componentList2) {
		final Set<ProteinComponent> intersection = getIntersection(componentList1, componentList2);
		final double numerator = 1.0 * intersection.size() * intersection.size();
		final double denominator = 1.0 * componentList1.size() * componentList2.size();
		return numerator / denominator;
	}

	public static double getSensitivity(Collection<ProteinComplex> predictedComplexes,
			Collection<ProteinComplex> referenceComplexes) {
		double numerator = 0.0;
		double denominator = 0.0;
		for (final ProteinComplex referenceComplex : referenceComplexes) {
			double t = 0.0;
			for (final ProteinComplex predictedComplex : predictedComplexes) {
				t = Math.max(t,
						Double.valueOf(getProteinComponentssInCommon(referenceComplex, predictedComplex).size()));
			}
			numerator += t;
			denominator += referenceComplex.getComponentList().size();
		}
		return numerator / denominator;
	}

	private static Set<ProteinComponent> getProteinComponentssInCommon(ProteinComplex complex1,
			ProteinComplex complex2) {

		final Set<ProteinComponent> intersection = new THashSet<ProteinComponent>();
		ProteinComplex smallComplex = complex1;
		ProteinComplex bigComplex = complex2;
		if (smallComplex.getComponentList().size() > bigComplex.getComponentList().size()) {
			smallComplex = complex2;
			bigComplex = complex1;
		}
		for (final ProteinComponent c1 : bigComplex.getComponentList()) {
			if (smallComplex.getComponentSet().contains(c1)) {
				intersection.add(c1);
			}
		}
		return intersection;
	}

	public static double getPositivePredictiveValue(Collection<ProteinComplex> predictedComplexes,
			Collection<ProteinComplex> referenceComplexes) {
		double numerator = 0.0;

		for (final ProteinComplex predictedComplex : predictedComplexes) {
			double t = 0.0;
			for (final ProteinComplex referenceComplex : referenceComplexes) {
				t = Math.max(t,
						Double.valueOf(getProteinComponentssInCommon(referenceComplex, predictedComplex).size()));
			}
			numerator += t;
		}
		double denominator = 0.0;
		for (final ProteinComplex predictedComplex : predictedComplexes) {
			for (final ProteinComplex referenceComplex : referenceComplexes) {
				denominator += getProteinComponentssInCommon(referenceComplex, predictedComplex).size();
			}
		}
		return numerator / denominator;
	}

	public static double getAccuracy(Collection<ProteinComplex> predictedComplexes,
			Collection<ProteinComplex> referenceComplexes) {
		final double ret = Math.sqrt(getSensitivity(predictedComplexes, referenceComplexes)
				* getPositivePredictiveValue(predictedComplexes, referenceComplexes));
		return ret;
	}

	public static double getMaximumMatchingRatio(Collection<ProteinComplex> predictedComplexes,
			Collection<ProteinComplex> referenceComplexes) {
		log.info("Calculating Maximum Matching Ratio between " + predictedComplexes.size() + " predicted complexes and "
				+ referenceComplexes.size() + " reference complexes");
		Map<ProteinComplex, Set<Edge>> edgesByReference = new THashMap<ProteinComplex, Set<Edge>>();
		Map<ProteinComplex, Set<Edge>> edgesByPredicted = new THashMap<ProteinComplex, Set<Edge>>();
		final List<Edge> edges = new ArrayList<Edge>();
		for (final ProteinComplex referenceComplex : referenceComplexes) {
			for (final ProteinComplex predictedComplex : predictedComplexes) {
				final double overlap = getOverlap(referenceComplex, predictedComplex);
				if (overlap > 0.0) {
					final Edge edge = new Edge(referenceComplex, predictedComplex, overlap);
					edges.add(edge);
					if (!edgesByReference.containsKey(referenceComplex)) {
						final Set<Edge> set = new THashSet<Edge>();
						edgesByReference.put(referenceComplex, set);
					}
					edgesByReference.get(referenceComplex).add(edge);
					if (!edgesByPredicted.containsKey(predictedComplex)) {
						final Set<Edge> set = new THashSet<Edge>();
						edgesByPredicted.put(predictedComplex, set);
					}
					edgesByPredicted.get(predictedComplex).add(edge);
				}
			}
		}
		// now we have to take the subset of edges so that their sum is maximum
		// and they each complex is not connected twice
		// therefore, the maximum number of elements of edges is the minimum
		// between the number of predicted and reference complexes that contain
		// edges
		final TDoubleList mmr = new TDoubleArrayList();
		filterEdges(edges, edgesByReference, edgesByPredicted);
		edgesByReference = getEdgesByReference(edges);
		edgesByPredicted = getEdgesByPredicted(edges);
		// first we divide the set of edges in sets of connected components
		final List<Set<Edge>> connectedComponents = getConnectedComponents(edges, edgesByPredicted, edgesByReference);
		for (final Set<Edge> listOfEdges : connectedComponents) {
			final Map<ProteinComplex, Set<Edge>> edgesByReference2 = getEdgesByReference(listOfEdges);
			final Map<ProteinComplex, Set<Edge>> edgesByPredicted2 = getEdgesByPredicted(listOfEdges);
			final int max = Math.min(edgesByPredicted2.size(), edgesByReference2.size());
			log.info("Maximum number of edges in the optimal set will be " + max);
			// now get all the combinations from 1 to max elements
			Double maxSum = -Double.MAX_VALUE;
			List<Edge> optimalEdgeCollection = null;
			for (int size = 1; size <= max; size++) {
				boolean someMaxFound = false;
				log.info("Getting possible collections of " + size + " size from " + listOfEdges.size() + " elements");
				final EdgeCombinationIterator edgeCollections = getCollectionsOfNElements(size, listOfEdges);
				// discard collections that connect to redundant complexes
				while (edgeCollections.hasNext()) {
					final List<Edge> edgeCollection = edgeCollections.next();
					final String progress = edgeCollections.getProgress();
					if (!"".equals(progress)) {
						log.info(progress);
					}
					if (edgeCollection == null) {
						continue;
					}
					if (!containRedundantComplexes(edgeCollection)) {
						final TDoubleList weights = new TDoubleArrayList();
						// get the sum of the weights of the edges
						edgeCollection.forEach(edge -> weights.add(edge.getWeight()));
						final double sumWeights = weights.sum();
						if (maxSum < sumWeights) {
							maxSum = sumWeights;
							optimalEdgeCollection = edgeCollection;
							someMaxFound = true;
						}
					} else {
						throw new IllegalArgumentException("This cannot happen");
					}
				}
				if (!someMaxFound) {
					break;
				}
			}
			// final double mmrTMP = maxSum / (optimalEdgeCollection.size() *
			// 1.0);
			final TDoubleList weights = new TDoubleArrayList();
			optimalEdgeCollection.forEach(edge -> weights.add(edge.getWeight()));
			mmr.addAll(weights);
		}

		final double mean = Maths.mean(mmr);
		log.info("MMR = " + mean);
		// now that we have the optimal weight collection, we calculate the MMR:
		return mean;
	}

	/**
	 * Discard edges when the complexes connected contain other edges with higher
	 * weights
	 * 
	 * @param edges
	 * @param edgesByPredicted
	 * @param edgesByReference
	 * @return
	 */
	private static void filterEdges(List<Edge> edges, Map<ProteinComplex, Set<Edge>> edgesByReference,
			Map<ProteinComplex, Set<Edge>> edgesByPredicted) {
		log.info("Filtering " + edges.size() + " edges to remove non informative ones");
		final Iterator<Edge> iterator = edges.iterator();
		int numRemoved = 0;
		while (iterator.hasNext()) {
			boolean remove1 = false;
			boolean remove2 = false;
			final Edge edge = iterator.next();
			final double w = edge.getWeight();
			final Set<Edge> edgesFromPredicted = edgesByPredicted.get(edge.getPredicted());
			for (final Edge edgeFromPredicted : edgesFromPredicted) {
				if (edgeFromPredicted != edge && edgeFromPredicted.getWeight() > w) {
					remove1 = true;
					break;
				}
			}

			final Set<Edge> edgesFromReference = edgesByReference.get(edge.getReference());
			for (final Edge edgeFromReference : edgesFromReference) {
				if (edgeFromReference != edge && edgeFromReference.getWeight() > w) {
					remove2 = true;
					break;
				}
			}
			if (remove1 && remove2) {
				iterator.remove();
				numRemoved++;
				continue;
			}
		}
		log.info("Number of edges removed: " + numRemoved);
	}

	private static List<Set<Edge>> getConnectedComponents(List<Edge> edges,
			Map<ProteinComplex, Set<Edge>> edgesByPredicted, Map<ProteinComplex, Set<Edge>> edgesByReference) {
		final Map<ProteinComplex, Set<Edge>> map = new THashMap<ProteinComplex, Set<Edge>>();
		for (final Edge edge : edges) {
			final ProteinComplex predicted = edge.getPredicted();
			if (!map.containsKey(predicted)) {
				map.put(predicted, new THashSet<Edge>());
			}
			if (edgesByPredicted.containsKey(predicted)) {
				map.get(predicted).addAll(edgesByPredicted.get(predicted));
			}

			final ProteinComplex reference = edge.getReference();
			if (!map.containsKey(reference)) {
				map.put(reference, new THashSet<Edge>());
			}
			if (edgesByReference.containsKey(reference)) {
				map.get(reference).addAll(edgesByReference.get(reference));
			}
		}

		final List<Set<Edge>> ret = new ArrayList<Set<Edge>>();
		for (final Set<Edge> list : map.values()) {
			ret.add(list);
		}

		log.info(ret.size());
		// now merge lists that have some element in common until no more
		// removal is done
		while (true) {
			boolean removal = false;
			for (int i = 0; i < ret.size(); i++) {
				final Set<Edge> list1 = ret.get(i);
				for (int j = i + 1; j < ret.size(); j++) {
					final Set<Edge> list2 = ret.get(j);
					if (containsElementsInCommon(list1, list2)) {
						list1.addAll(list2);
						list2.clear();
					}
				}
			}
			log.info(ret.size());
			// remove empty sets
			final Iterator<Set<Edge>> iterator = ret.iterator();
			while (iterator.hasNext()) {
				final Set<Edge> next = iterator.next();
				if (next.isEmpty()) {
					iterator.remove();
					removal = true;
				}
			}
			if (!removal) {
				break;
			}
		}
		log.info(ret.size());
		return ret;
	}

	private static boolean containsElementsInCommon(Set<Edge> list1, Set<Edge> list2) {
		for (final Edge edge : list2) {
			if (list1.contains(edge)) {
				return true;
			}
		}
		return false;
	}

	private static Map<ProteinComplex, Set<Edge>> getEdgesByReference(Collection<Edge> listOfEdges) {
		final Map<ProteinComplex, Set<Edge>> edgesByReference = new THashMap<ProteinComplex, Set<Edge>>();
		for (final Edge edge : listOfEdges) {
			final ProteinComplex reference = edge.getReference();
			if (!edgesByReference.containsKey(reference)) {
				edgesByReference.put(reference, new THashSet<Edge>());
			}
			edgesByReference.get(reference).add(edge);
		}
		return edgesByReference;
	}

	private static Map<ProteinComplex, Set<Edge>> getEdgesByPredicted(Collection<Edge> listOfEdges) {
		final Map<ProteinComplex, Set<Edge>> edgesByPredicted = new THashMap<ProteinComplex, Set<Edge>>();
		for (final Edge edge : listOfEdges) {
			final ProteinComplex predicted = edge.getPredicted();
			if (!edgesByPredicted.containsKey(predicted)) {
				edgesByPredicted.put(predicted, new THashSet<Edge>());
			}
			edgesByPredicted.get(predicted).add(edge);
		}
		return edgesByPredicted;
	}

	/**
	 * Returns all the combinations of N elements from the list of objects
	 * 
	 * @param n
	 * @param objs
	 * @return
	 */
	private static EdgeCombinationIterator getCollectionsOfNElements(int size, Collection<Edge> edges) {

		return new EdgeCombinationIterator(size, edges);

	}

	private static boolean containRedundantComplexes(List<Edge> edgeCollection) {
		final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
		for (final Edge edge : edgeCollection) {
			final ProteinComplex predicted = edge.getPredicted();
			if (set.contains(predicted)) {
				return true;
			}
			set.add(predicted);
			final ProteinComplex reference = edge.getReference();
			if (set.contains(reference)) {
				return true;
			}
			set.add(reference);
		}
		return false;
	}

	public static double calculateCompositeScore(double maxMatchingRatio, List<ProteinComplex> predictedComplexes,
			Set<ProteinComplex> referenceComplexes) {

		final double accuracy = getAccuracy(predictedComplexes, referenceComplexes);
		final double compositeScore = maxMatchingRatio + accuracy;
		return compositeScore;
	}
}
