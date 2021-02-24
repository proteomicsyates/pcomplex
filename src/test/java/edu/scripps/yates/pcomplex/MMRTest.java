package edu.scripps.yates.pcomplex;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import gnu.trove.map.hash.THashMap;
import junit.framework.Assert;

public class MMRTest {
	private final static Map<String, ProteinComponent> proteins = new THashMap<String, ProteinComponent>();

	@Test
	public void mmrCalculationTest() throws IOException {
		List<ProteinComplex> referenceComplexes = createReference(false);
		List<ProteinComplex> predictedComplexes = createPredictedComplexes(true);
		double maximumMatchingRatio = ClusterEvaluation.getMaximumMatchingRatio(predictedComplexes, referenceComplexes);
		Assert.assertEquals(0.775, maximumMatchingRatio);

		referenceComplexes = createReference(true);
		predictedComplexes = createPredictedComplexes(true);
		maximumMatchingRatio = ClusterEvaluation.getMaximumMatchingRatio(predictedComplexes, referenceComplexes);
		Assert.assertEquals(0.738888888888, maximumMatchingRatio, 0.000001);
	}

	private List<ProteinComplex> createPredictedComplexes(boolean b) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final String[] p1proteins = { "b", "c", "d", "e" };
		final ProteinComplex p1 = createComplex("P1", p1proteins);
		ret.add(p1);
		final String[] p2proteins = { "a", "b", "c", "f" };
		final ProteinComplex p2 = createComplex("P2", p2proteins);
		ret.add(p2);
		final String[] p3proteins = { "f", "g", "h" };
		final ProteinComplex p3 = createComplex("P3", p3proteins);
		ret.add(p3);
		if (b) {
			final String[] p4proteins = { "j", "k", "l" };
			final ProteinComplex p4 = createComplex("P4", p4proteins);
			ret.add(p4);
		}
		return ret;
	}

	private List<ProteinComplex> createReference(boolean b) throws IOException {
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final String[] r1proteins = { "a", "b", "c", "d", "e" };
		final ProteinComplex r1 = createComplex("R1", r1proteins);
		ret.add(r1);
		final String[] r2proteins = { "f", "g", "h", "i" };
		final ProteinComplex r2 = createComplex("R2", r2proteins);
		ret.add(r2);
		if (b) {
			final String[] r3proteins = { "j", "l" };
			final ProteinComplex r3 = createComplex("R3", r3proteins);
			ret.add(r3);
		}
		return ret;
	}

	private ProteinComplex createComplex(String id, String[] proteinNames) throws IOException {
		final ProteinComplex ret = new ProteinComplex(id);
		for (final String proteinName : proteinNames) {
			ret.addComponent(getProteinComponent(proteinName));
		}
		return ret;
	}

	private ProteinComponent getProteinComponent(String proteinName) throws IOException {
		if (proteins.containsKey(proteinName)) {
			return proteins.get(proteinName);
		}
		final ProteinComponent component = new ProteinComponent(proteinName, null);
		proteins.put(proteinName, component);
		return component;
	}
}
