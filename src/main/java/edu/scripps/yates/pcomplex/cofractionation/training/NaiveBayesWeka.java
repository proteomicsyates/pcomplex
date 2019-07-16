package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.DistanceMeasure;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import smile.netlib.NLMatrix;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.NaiveBayes;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffLoader;

public abstract class NaiveBayesWeka extends NaiveBayesClassification {
	private final static Logger log = Logger.getLogger(NaiveBayesWeka.class);
	private final String datasetName;
	private final String arfFileFullPath;

	public NaiveBayesWeka(String datasetName, String arfFileFullPath) {
		super();
		this.datasetName = datasetName;
		this.arfFileFullPath = arfFileFullPath;
	}

	@Override
	public abstract ClassLabel getClassLabel(String protein, String protein2) throws IOException;

	@Override
	public ClassificationResult naiveBayesTraining(Map<DistanceMeasure, NLMatrix> matrixMap, List<String> proteinList)
			throws Exception {
		log.info("Starting Naive Bayesian approach");
		final File arffFile = createArffFile(datasetName, matrixMap, proteinList);
		final ArffLoader loader = new ArffLoader();
		loader.setSource(arffFile);
		final Instances dataSet = loader.getDataSet();
		final int classIndex = dataSet.numAttributes() - 1;
		dataSet.setClassIndex(classIndex);
		final Classifier classifier = new NaiveBayes();
		final int kFold = 15;
		for (int foldNumber = 0; foldNumber < kFold; foldNumber++) {
			final Instances trainCV = dataSet.trainCV(kFold, foldNumber);

			classifier.buildClassifier(trainCV);

			final Instances testCV = dataSet.testCV(kFold, foldNumber);
			final Evaluation eval = new Evaluation(trainCV);
			// final AbstractOutput forPredictionsPrinting = new
			// AbstractOutput() {
			//
			// @Override
			// public String globalInfo() {
			// // TODO Auto-generated method stub
			// return null;
			// }
			//
			// @Override
			// public String getDisplay() {
			// // TODO Auto-generated method stub
			// return null;
			// }
			//
			// @Override
			// protected void doPrintHeader() {
			// // TODO Auto-generated method stub
			//
			// }
			//
			// @Override
			// protected void doPrintFooter() {
			// // TODO Auto-generated method stub
			//
			// }
			//
			// @Override
			// protected void doPrintClassification(double[] dist, Instance
			// inst, int index) throws Exception {
			// // TODO Auto-generated method stub
			//
			// }
			//
			// @Override
			// protected void doPrintClassification(Classifier classifier,
			// Instance inst, int index) throws Exception {
			// // TODO Auto-generated method stub
			//
			// }
			// };
			eval.evaluateModel(classifier, testCV);// , forPredictionsPrinting);
			System.out.println("** Naive Bayes Evaluation with Datasets **");
			System.out.println(eval.toSummaryString());
			System.out.print(" the expression for the input data as per alogorithm is ");
			System.out.println(classifier);

		}
		final ClassificationResult result = new ClassificationResult(proteinList.size(), proteinList.size(),
				ClassLabel.values().length);
		int instanceIndex = 0;
		for (int i = 0; i < proteinList.size(); i++) {
			for (int j = i + 1; j < proteinList.size(); j++) {
				final ClassLabel classLabel = getClassLabel(proteinList.get(i), proteinList.get(j));
				final Instance instance = dataSet.get(instanceIndex);
				final String value = instance.attribute(0).value(0);
				final double[] probabilities = classifier.distributionForInstance(instance);
				result.set(i, j, probabilities, classLabel);
				instanceIndex++;
			}
		}
		return result;
	}

	/**
	 * Creates a tmp file in arff format, following description on
	 * https://machinelearningmastery.com/load-csv-machine-learning-data-weka/
	 * 
	 * @param matrixMap
	 * @param proteinList
	 * @return
	 * @throws IOException
	 */
	private File createArffFile(String datasetName, Map<DistanceMeasure, NLMatrix> matrixMap, List<String> proteinList)
			throws IOException {
		File file = null;
		if (arfFileFullPath != null) {
			file = new File(arfFileFullPath);
		} else {
			file = File.createTempFile("measurements", "arff");
			file.deleteOnExit();
		}
		final FileWriter fw = new FileWriter(file);
		// write header
		fw.write("@RELATION " + datasetName + "\n\n");
		for (final DistanceMeasure distance : DistanceMeasure.values()) {
			fw.write("@ATTRIBUTE " + distance.name() + " REAL\n");
		}
		fw.write("@ATTRIBUTE class {");
		final StringBuilder sb2 = new StringBuilder();
		for (final ClassLabel classLabel : ClassLabel.values()) {
			if (!"".equals(sb2.toString())) {
				sb2.append(",");
			}
			sb2.append(classLabel.name());
		}
		fw.write(sb2.toString() + "}\n\n");
		fw.write("@DATA\n");
		// write data
		for (int i = 0; i < proteinList.size(); i++) {
			for (int j = i + 1; j < proteinList.size(); j++) {
				final ClassLabel classLabel = getClassLabel(proteinList.get(i), proteinList.get(j));

				final StringBuilder sb = new StringBuilder();
				for (final DistanceMeasure distance : DistanceMeasure.values()) {
					final NLMatrix nlMatrix = matrixMap.get(distance);
					final double num = nlMatrix.get(i, j);
					if (!"".equals(sb.toString())) {
						sb.append(",");
					}
					sb.append(num);
				}
				sb.append(",").append(classLabel.name());
				fw.write(sb.toString() + "\n");
			}
		}
		fw.close();
		return file;
	}

}
