package edu.scripps.yates.pcomplex.util;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.Enzyme;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

public class StochiometryEstimator {

	private final UniprotProteinLocalRetriever uplr;
	private final List<String> accs = new ArrayList<String>();
	private final SeparationExperiment exp;
	private boolean loaded = false;
	private final Enzyme enzyme;
	private static int MIN_FRACTIONS = 10;

	/**
	 * 
	 * @param uniprotFolder
	 *            to get the protein sequences
	 * @param exp
	 * @param enzyme
	 *            to make the in-silico digestion
	 */
	public StochiometryEstimator(File uniprotFolder, SeparationExperiment exp, Enzyme enzyme) {
		uplr = new UniprotProteinLocalRetriever(uniprotFolder, true);
		accs.addAll(exp.getProteinACCList());
		this.exp = exp;
		this.enzyme = enzyme;
	}

	private void loadAnnotations() {
		if (!loaded) {
			uplr.getAnnotatedProteins(null, accs);
			loaded = true;
		}
	}

	public double getStochiometry(String protein1, String protein2) {
		if (!loaded) {
			loadAnnotations();
		}
		// only can be estimated if both proteins are present in at least
		// MIN_FRACTIONS
		final TDoubleArrayList profile1 = exp.getElutionProfile(protein1, DataType.SPC);
		int numPresent = 0;
		for (final double spc : profile1.toArray()) {
			if (!Double.isNaN(spc) && spc > 0.0) {
				numPresent++;
			}
		}
		if (numPresent < MIN_FRACTIONS) {
			return Double.NaN;
		}
		final TDoubleArrayList profile2 = exp.getElutionProfile(protein2, DataType.SPC);
		numPresent = 0;
		for (final double spc : profile2.toArray()) {
			if (!Double.isNaN(spc) && spc > 0.0) {
				numPresent++;
			}
		}
		if (numPresent < MIN_FRACTIONS) {
			return Double.NaN;
		}
		// now, calculate the estimated number of peptides
		final int e1 = getNumTheoreticalSPC(protein1);
		final int e2 = getNumTheoreticalSPC(protein2);
		if (e1 > 0.0 && e2 > 0.0) {
			final TDoubleList nums = new TDoubleArrayList();
			for (int i = 0; i < profile1.size(); i++) {
				final double spc1 = profile1.get(i);
				final double spc2 = profile2.get(i);
				if (!Double.isNaN(spc1) && !Double.isNaN(spc2) && spc1 > 0.0 && spc2 > 0.0) {
					final double numerator = spc1 / e1;
					final double denominator = spc1 / e1;
					final double ratio = numerator / denominator;
					nums.add(ratio);
				}
			}
			final Median median = new Median();
			final double stochiometry = median.evaluate(nums.toArray());
			return stochiometry;
		}
		return Double.NaN;
	}

	private int getNumTheoreticalSPC(String protein1) {
		final Map<String, Entry> entryMap = uplr.getAnnotatedProtein(null, protein1);
		if (entryMap.containsKey(protein1)) {
			final String sequence = UniprotEntryUtil.getProteinSequence(entryMap.get(protein1));
			final Enzyme trypsin = getEnzyme();

			final Protein protein = new Protein(new AASequenceImpl(sequence));
			final Protein[] cleave = trypsin.cleave(protein);
			return cleave.length;
		}
		return -1;
	}

	private Enzyme getEnzyme() {
		return enzyme;
	}
}
