package edu.scripps.yates.pcomplex;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.cofractionation.training.ClassificationResult;
import edu.scripps.yates.pcomplex.cofractionation.training.ProteinPairInteraction;
import edu.scripps.yates.pcomplex.model.Cluster;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.utilities.dates.DatesUtil;
import gnu.trove.map.hash.TObjectIntHashMap;
import uk.ac.rhul.cs.cl1.ClusterONE;
import uk.ac.rhul.cs.cl1.ClusterONEAlgorithmParameters;
import uk.ac.rhul.cs.cl1.ClusterONEException;
import uk.ac.rhul.cs.cl1.TaskMonitor;
import uk.ac.rhul.cs.cl1.ValuedNodeSet;
import uk.ac.rhul.cs.graph.Graph;

public class ClusterOneInterface {
	private final static Logger log = Logger.getLogger(Class.class);
	private final TObjectIntHashMap<String> nodesByProtein = new TObjectIntHashMap<String>();

	public List<ProteinComplex> runClusterOne(ClassificationResult classificationResult, double pValueCutoff,
			ClassLabel classLabel, double minDensity, double penalty) {
		final ClusterONEAlgorithmParameters algorithmParameters = new ClusterONEAlgorithmParameters();
		algorithmParameters.setMinDensity(minDensity);
		algorithmParameters.setNodePenalty(penalty);

		final ClusterONE clusterOne = new ClusterONE(algorithmParameters);
		clusterOne.setGraph(createGraph(classificationResult, pValueCutoff, classLabel));
		clusterOne.setTaskMonitor(new TaskMonitor() {

			@Override
			public void setStatus(String message) {
				log.info("ClusterONE status: " + message);
			}

			@Override
			public void setPercentCompleted(int percent) throws IllegalArgumentException {
				log.info("ClusterONE progress: " + percent);

			}

			@Override
			public void setException(Throwable t, String userErrorMessage, String recoveryTip) {
				log.error(userErrorMessage, t);
				log.error("ClusterONE tip: " + recoveryTip);
			}

			@Override
			public void setException(Throwable t, String userErrorMessage) {
				log.error(userErrorMessage, t);
			}

			@Override
			public void setEstimatedTimeRemaining(long time) {
				log.info("ClusterONE ETA: " + DatesUtil.getDescriptiveTimeFromMillisecs(time));
			}
		});
		try {
			clusterOne.run();
			final List<ValuedNodeSet> results = clusterOne.getResults();

			return trasnformResultsToProteinComplexes(results);
		} catch (final ClusterONEException e) {

			e.printStackTrace();
			log.error(e);
			return null;
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e);
			return null;
		}
	}

	private List<ProteinComplex> trasnformResultsToProteinComplexes(List<ValuedNodeSet> results) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		for (final ValuedNodeSet nodeSet : results) {
			ret.add(createProteinComplex(nodeSet));
		}
		return ret;
	}

	private ProteinComplex createProteinComplex(ValuedNodeSet nodeSet) throws IOException {

		final Cluster cluster = new Cluster(nodeSet);
		final ProteinComplex complex = new ProteinComplex(String.valueOf(cluster.hashCode()));
		final List<ProteinComponent> proteinComponents = cluster.getProteinComponents();
		for (final ProteinComponent proteinComponent : proteinComponents) {
			complex.addComponent(proteinComponent);
			complex.setDensity(cluster.getDensity());
			complex.setSignificance(cluster.getSignificance());
		}
		return complex;
	}

	/**
	 * Creates the {@link Graph}
	 * 
	 * @param result
	 * @param pValueCutoff
	 * @return
	 */
	private Graph createGraph(ClassificationResult result, double pValueCutoff, ClassLabel classLabel) {
		final Graph ret = new Graph();
		final List<ProteinPairInteraction> interactionsByPValue = result.getInteractions(pValueCutoff, classLabel);
		int edges = 0;
		for (final ProteinPairInteraction interaction : interactionsByPValue) {
			if (interaction.getProbability(ClassLabel.INTRA_COMPLEX) < pValueCutoff) {
				throw new IllegalArgumentException(
						"This shoudn't happen because it should be already filtered by getInteractionsByPValue()");
			}
			final String protein1 = interaction.getProtein1Acc();
			final String protein2 = interaction.getProtein2Acc();
			int indexNode1 = -1;
			if (nodesByProtein.containsKey(protein1)) {
				indexNode1 = nodesByProtein.get(protein1);
			} else {
				indexNode1 = ret.createNode(protein1);
				nodesByProtein.put(protein1, indexNode1);
			}
			int indexNode2 = -1;
			if (nodesByProtein.containsKey(protein2)) {
				indexNode2 = nodesByProtein.get(protein2);
			} else {
				indexNode2 = ret.createNode(protein2);
				nodesByProtein.put(protein2, indexNode2);
			}
			if (ret.areConnected(indexNode1, indexNode2)) {
				log.info("This shouldnt happen?");
			} else {
				final double weight = interaction.getProbability(ClassLabel.INTRA_COMPLEX);
				ret.createEdge(indexNode1, indexNode2, weight);
				edges++;
			}
		}
		log.info("Graph created with " + edges + " edges");
		if (edges == 0) {
			throw new IllegalArgumentException(
					"There is no edges in this graph formed with " + interactionsByPValue.size() + " interactions");
		}
		return ret;
	}

}
