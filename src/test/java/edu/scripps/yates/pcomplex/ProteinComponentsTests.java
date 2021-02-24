package edu.scripps.yates.pcomplex;

import java.io.IOException;
import java.util.Set;

import org.junit.Test;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;
import gnu.trove.set.hash.THashSet;
import junit.framework.Assert;

public class ProteinComponentsTests {
	@Test
	public void proteinComponentComparisonTest() throws IOException {
		final String A = "A";
		final ProteinComponent c1 = new ProteinComponent(A, A);
		final String B = "B" + ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR + "A";
		final ProteinComponent c2 = new ProteinComponent(B, B);
		Assert.assertEquals(c1, c2);
	}

	@Test
	public void proteinComplexComparisonTest() throws IOException {
		final String A = "A";
		final ProteinComponent c1 = new ProteinComponent(A, A);
		final String B = "B" + ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR + "A";
		final ProteinComponent c2 = new ProteinComponent(B, B);
		final ProteinComplex complex1 = new ProteinComplex("a");
		complex1.addComponent(c1);

		final ProteinComplex complex2 = new ProteinComplex("b");
		complex2.addComponent(c2);

		Assert.assertEquals(complex1, complex2);
	}

	@Test
	public void proteinComplexOverlapTest() throws IOException {
		final String A = "A";
		final ProteinComponent c1 = new ProteinComponent(A, A);
		final String B = "B" + ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR + "A";
		final ProteinComponent c2 = new ProteinComponent(B, B);
		final ProteinComplex complex1 = new ProteinComplex("a");
		complex1.addComponent(c1);

		final ProteinComplex complex2 = new ProteinComplex("b");
		complex2.addComponent(c2);

		final double overlap = ClusterEvaluation.getOverlap(complex1, complex2);
		Assert.assertEquals(1.0, overlap);
	}

	@Test
	public void proteinComplexSetOperations() throws IOException {
		final Set<ProteinComplex> set = new THashSet<ProteinComplex>();
		final String A = "A";
		final ProteinComponent c1 = new ProteinComponent(A, A);
		final String B = "A" + ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR + "B";
		final ProteinComponent c2 = new ProteinComponent(B, B);

		final ProteinComplex complex1 = new ProteinComplex("a");
		complex1.addComponent(c1);

		final ProteinComplex complex2 = new ProteinComplex("b");
		complex2.addComponent(c2);

		set.add(complex1);
		set.add(complex2);
		Assert.assertEquals(1, set.size()); // only 1 because they are the same

		final ProteinComplex complex3 = new ProteinComplex("c");
		final String C = "C";
		final ProteinComponent c3 = new ProteinComponent(C, C);
		complex3.addComponent(c3);
		set.add(complex3);
		Assert.assertEquals(2, set.size());

		// remove c1
		set.remove(complex1);
		Assert.assertEquals(1, set.size());
		Assert.assertFalse(set.contains(complex1));
		Assert.assertFalse(set.contains(complex2));
		Assert.assertTrue(set.contains(complex3));
	}
}
