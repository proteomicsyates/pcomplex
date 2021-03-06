package edu.scripps.yates.pcomplex;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import junit.framework.Assert;
import smile.classification.FLD;
import smile.data.AttributeDataset;
import smile.data.parser.ArffParser;
import smile.validation.LOOCV;

public class tests {
	@Test
	public void test1() throws IOException {
		final ProteinComplex complex1 = new ProteinComplex("ASDF");
		complex1.addComponent(new ProteinComponent("COMP1", ""));
		complex1.addComponent(new ProteinComponent("COMP2", ""));
		System.out.println(complex1.hashCode());

		final ProteinComplex complex2 = new ProteinComplex("ASDF12");
		complex2.addComponent(new ProteinComponent("COMP2", ""));
		complex2.addComponent(new ProteinComponent("COMP1", ""));
		System.out.println(complex2.hashCode());

		Assert.assertEquals(complex1.hashCode(), complex2.hashCode());
		complex2.addComponent(new ProteinComponent("COMP3", ""));
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

	@Test
	public void dtaselecttest() {
		DTASelectParser parser;
		try {
			parser = new DTASelectParser(new File(
					"C:\\Users\\salvador\\Desktop\\Anthony\\protein_complexes\\experiments\\Beta_cell_PCP\\106-PBC_106_S3_A9_1_3389_nopd_2020_10_28_08_255851-174397-DTASelect.txt"));
			parser.getPSMsByFullSequence();
			final String fastaPath = parser.getFastaPath();
			System.out.println(fastaPath);
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}
}
