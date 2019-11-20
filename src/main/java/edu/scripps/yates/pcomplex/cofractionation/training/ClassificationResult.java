package edu.scripps.yates.pcomplex.cofractionation.training;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import gnu.trove.map.hash.THashMap;

public class ClassificationResult {

	private final Map<String, Map<String, ProteinPairInteraction>> interactionsByProteinAcc = new THashMap<String, Map<String, ProteinPairInteraction>>();
	private final List<ProteinPairInteraction> proteinPairInteractions = new ArrayList<ProteinPairInteraction>();
	private final boolean precisionsAreCalculated = false;
	private final MyClassifier myClassifier;

	public ClassificationResult(MyClassifier myClassifier) {
		this.myClassifier = myClassifier;
	}

	public void addProteinInteraction(int i, int j, List<String> proteinList, double[] probabilities,
			ClassLabel trueLabel) {
		final ProteinPairInteraction ppi = new ProteinPairInteraction(proteinList.get(i), proteinList.get(j),
				probabilities, trueLabel);
		addProteinInteraction(ppi);
	}

	public void addProteinInteraction(ProteinPairInteraction ppi) {

		proteinPairInteractions.add(ppi);
		// add by protein1
		if (!interactionsByProteinAcc.containsKey(ppi.getProtein1Acc())) {
			interactionsByProteinAcc.put(ppi.getProtein1Acc(), new THashMap<String, ProteinPairInteraction>());
		}
		interactionsByProteinAcc.get(ppi.getProtein1Acc()).put(ppi.getProtein2Acc(), ppi);
		// add by protein2
		if (!interactionsByProteinAcc.containsKey(ppi.getProtein2Acc())) {
			interactionsByProteinAcc.put(ppi.getProtein2Acc(), new THashMap<String, ProteinPairInteraction>());
		}
		interactionsByProteinAcc.get(ppi.getProtein2Acc()).put(ppi.getProtein1Acc(), ppi);
	}

	/**
	 * Gets the {@link ProteinPairInteraction} with a pValue >= the cutoff
	 * provided for the provided {@link ClassLabel}
	 * 
	 * @param pValueCutOff
	 * @param intraComplex
	 * @return
	 */
	public List<ProteinPairInteraction> getInteractions(double pValueCutOff, ClassLabel classLabel) {
		final List<ProteinPairInteraction> ret = new ArrayList<ProteinPairInteraction>();
		for (final ProteinPairInteraction ppi : proteinPairInteractions) {
			final double probability = ppi.getProbability(classLabel);
			if (probability >= pValueCutOff) {
				ret.add(ppi);
			}
		}
		return ret;
	}

	public MyClassifier getMyClassifier() {
		return myClassifier;
	}

}
