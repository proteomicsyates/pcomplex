package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.IOException;
import java.util.Collection;

import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;

public interface TrueClassifier {

	public ClassLabel getClassLabel(String protein1, String protein2) throws IOException;

	public ClassLabel getClassLabel(Collection<String> proteins, Collection<String> proteins2) throws IOException;

}
