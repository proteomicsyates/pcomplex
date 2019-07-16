package edu.scripps.yates.pcomplex;

import static org.junit.Assert.assertEquals;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import junit.framework.Assert;
import smile.classification.FLD;
import smile.data.AttributeDataset;
import smile.data.parser.ArffParser;
import smile.validation.LOOCV;

public class tests {
	@Test
	public void test1() {
		final ProteinComplex complex1 = new ProteinComplex("ASDF", true);
		complex1.addComponent("COMP1");
		complex1.addComponent("COMP2");
		System.out.println(complex1.hashCode());

		final ProteinComplex complex2 = new ProteinComplex("ASDF12", true);
		complex2.addComponent("COMP2");
		complex2.addComponent("COMP1");
		System.out.println(complex2.hashCode());

		Assert.assertEquals(complex1.hashCode(), complex2.hashCode());
		complex2.addComponent("COMP3");
		System.out.println(complex2.hashCode());
		Assert.assertNotSame(complex1.hashCode(), complex2.hashCode());

	}

	@Test
	public void fisher() {
		System.out.println("IRIS");
		final ArffParser arffParser = new ArffParser();
		arffParser.setResponseIndex(4);
		try {
			final AttributeDataset iris = arffParser.parse(new ClassPathResource("iris.arff").getFile());
			final double[][] x = iris.toArray(new double[iris.size()][]);
			final int[] y = iris.toArray(new int[iris.size()]);

			final int n = x.length;
			final LOOCV loocv = new LOOCV(n);
			int error = 0;
			for (int i = 0; i < n; i++) {
				final double[][] trainx = smile.math.Math.slice(x, loocv.train[i]);
				final int[] trainy = smile.math.Math.slice(y, loocv.train[i]);
				final FLD fisher = new FLD(trainx, trainy);

				if (y[loocv.test[i]] != fisher.predict(x[loocv.test[i]]))
					error++;
			}

			System.out.println("FLD error = " + error);
			assertEquals(5, error);
		} catch (final Exception ex) {
			System.err.println(ex);
		}
	}
}
