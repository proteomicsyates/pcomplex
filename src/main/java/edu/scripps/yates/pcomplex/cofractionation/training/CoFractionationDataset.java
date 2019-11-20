package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.DistanceMeasure;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import smile.netlib.NLMatrix;

public class CoFractionationDataset {
	private final static Logger log = Logger.getLogger(CoFractionationDataset.class);

	private final Dataset dataset;
	private final Dataset trainingDataset;

	public CoFractionationDataset(Map<DistanceMeasure, NLMatrix> distances, List<String> proteinKeys,
			TrueClassifier trueClassifier, double minCorrelationForTraining, ProteinComplexDB referenceDB, File folder,
			String datasetName) throws IOException {

		// no filters
		dataset = new Dataset(false, datasetName, folder, trueClassifier, distances, proteinKeys);
		// filtering by correlation and by referenceDB
		trainingDataset = new Dataset(true, datasetName, folder, trueClassifier, distances, proteinKeys,
				minCorrelationForTraining, referenceDB);
	}

	public Dataset getTrainingDataset() {
		return trainingDataset;
	}

	public Dataset getDataset() {
		return dataset;
	}

}
