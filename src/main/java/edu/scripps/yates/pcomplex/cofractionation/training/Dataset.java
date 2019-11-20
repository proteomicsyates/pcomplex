package edu.scripps.yates.pcomplex.cofractionation.training;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.DistanceMeasure;
import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.cofractionation.ClassLabel;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.THashSet;
import smile.netlib.NLMatrix;
import weka.core.Instances;
import weka.core.converters.ArffLoader;

public class Dataset {
	private final static Logger log = Logger.getLogger(Dataset.class);
	private final Map<DistanceMeasure, TObjectDoubleMap<String>> distances = new THashMap<DistanceMeasure, TObjectDoubleMap<String>>();
	private final List<String> proteins = new ArrayList<String>();
	private Instances instances;
	private final String name;
	private final File folder;
	private List<DistanceMeasure> distanceList;
	private final List<String> proteinPairKeysList;
	private static final String PAIR_SEPARATOR = "@@@@";
	private final TrueClassifier trueClassifier;
	private final boolean training;

	public Dataset(boolean training, String name, File folder, TrueClassifier trueClassifier,
			Map<DistanceMeasure, NLMatrix> distances, List<String> proteinKeys) {
		this(training, name, folder, trueClassifier, distances, proteinKeys, -Double.MAX_VALUE, null);
	}

	public Dataset(boolean training, String name, File folder, TrueClassifier trueClassifier,
			Map<DistanceMeasure, NLMatrix> distances, List<String> proteinKeys, double minCorrelation) {
		this(training, name, folder, trueClassifier, distances, proteinKeys, minCorrelation, null);
	}

	public Dataset(boolean training, String name, File folder, TrueClassifier trueClassifier,
			Map<DistanceMeasure, NLMatrix> distances, List<String> proteinKeys, ProteinComplexDB referenceDB) {
		this(training, name, folder, trueClassifier, distances, proteinKeys, -Double.MAX_VALUE, referenceDB);
	}

	public Dataset(boolean training, String name, File folder, TrueClassifier trueClassifier,
			Map<DistanceMeasure, NLMatrix> distances, List<String> proteinKeys, double minCorrelation,
			ProteinComplexDB referenceDB) {
		this.training = training;
		this.name = name;
		this.folder = folder;
		this.trueClassifier = trueClassifier;
		final NLMatrix correlations = distances.get(DistanceMeasure.PEARSON_CORRELATION_COEFFICIENT);

		final Set<String> proteinPairKeys = new THashSet<String>();
		proteinPairKeysList = new ArrayList<String>();
		for (int index = 0; index < proteinKeys.size(); index++) {
			final String protein1 = proteinKeys.get(index);
			if (referenceDB != null && referenceDB
					.getProteinComplexesByProteins(getIndividualProteinsFromProteinGroupKey(protein1)).isEmpty()) {
				continue;
			}
			String protein2 = null;
			boolean proteinPairValid = false;
			for (int index2 = 0; index2 < proteinKeys.size(); index2++) {
				protein2 = proteinKeys.get(index2);
				if (referenceDB != null && referenceDB
						.getProteinComplexesByProteins(getIndividualProteinsFromProteinGroupKey(protein2)).isEmpty()) {
					continue;
				}

				final double correlation = correlations.get(index, index2);
				if (correlation >= minCorrelation) {
					proteinPairValid = true;
					final String pairKey = getPairKey(protein1, protein2);
					if (!proteinPairKeys.contains(pairKey)) {
						proteinPairKeys.add(pairKey);
						proteinPairKeysList.add(pairKey);
					}
					for (final DistanceMeasure distance : distances.keySet()) {
						if (!this.distances.containsKey(distance)) {
							this.distances.put(distance, new TObjectDoubleHashMap<String>());
						}
						final TObjectDoubleMap<String> valueMap = this.distances.get(distance);
						final NLMatrix matrix = distances.get(distance);
						valueMap.put(pairKey, matrix.get(index, index2));
					}
				}
			}
			if (proteinPairValid) {
				proteins.add(protein1);
				proteins.add(protein2);
			}
		}
		log.info("Dataset " + name + " created with " + getNumProteins() + " proteins and " + proteinPairKeys.size()
				+ " protein pairs");
	}

	private List<String> getIndividualProteinsFromProteinGroupKey(String proteinGroupKey) {
		final List<String> ret = new ArrayList<String>();
		if (proteinGroupKey.contains(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR)) {
			final String[] split = proteinGroupKey.split(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
			for (final String acc : split) {
				ret.add(acc);
			}
		} else {
			ret.add(proteinGroupKey);
		}
		return ret;
	}

	public List<String> getProteinKeys() {
		return proteins;
	}

	public String getProteinKey(int index) {
		return proteins.get(index);
	}

	public int getNumProteins() {
		return proteins.size();
	}

	public double getValue(DistanceMeasure distance, String protein1, String protein2) {
		return getValue(distance, getPairKey(protein1, protein2));
	}

	public double getValue(DistanceMeasure distance, String proteinPairKey) {
		if (!distances.containsKey(distance)) {
			return Double.NaN;
		}
		final TObjectDoubleMap<String> valuesByProteinPairKey = distances.get(distance);
		final double d = valuesByProteinPairKey.get(proteinPairKey);
		return d;
	}

	private String getPairKey(String protein1, String protein2) {
		return protein1 + PAIR_SEPARATOR + protein2;
	}

	public List<DistanceMeasure> getDistances() {
		if (distanceList == null) {
			distanceList = new ArrayList<DistanceMeasure>();

			distanceList.addAll(distances.keySet());
			Collections.sort(distanceList, new Comparator<DistanceMeasure>() {

				@Override
				public int compare(DistanceMeasure o1, DistanceMeasure o2) {
					return Integer.compare(o1.ordinal(), o2.ordinal());
				}
			});
		}
		return distanceList;
	}

	public int getNumInstances() {
		// because we have everything repeated with symmetric keys we divided by
		// 2
		return distances.get(distances.keySet().iterator().next()).size() / 2;
	}

	public Instances getInstances() throws IOException {

		if (instances == null) {
			String pathname = folder.getAbsolutePath() + File.separator + name;
			if (isTraining()) {
				pathname += "_train";
			}
			pathname += ".arff";
			final File arffFile = new File(pathname);
			createArffFile(arffFile, trueClassifier, isTraining());
			final ArffLoader loader = new ArffLoader();
			loader.setSource(arffFile);
			instances = loader.getDataSet();
			if (isTraining()) {
				instances.setClassIndex(instances.numAttributes() - 1);
				final int numClasses = instances.numClasses();
				log.info(instances.size() + " interactions with " + numClasses + " classes in ARFF file "
						+ arffFile.getAbsolutePath());
			} else {
				log.info(instances.size() + " interactions in ARFF file " + arffFile.getAbsolutePath());
			}

		}
		return instances;
	}

	private boolean isTraining() {
		return training;
	}

	/**
	 * Creates a ARFF file following description on
	 * https://machinelearningmastery.com/load-csv-machine-learning-data-weka/
	 * 
	 * @param arffFile
	 * @param dataset
	 * @param minCorrelation
	 * @param trueClassifier
	 * @param writeClass
	 * @throws IOException
	 */
	private void createArffFile(File arffFile, TrueClassifier trueClassifier, boolean writeClass) throws IOException {

		final FileWriter fw = new FileWriter(arffFile);
		// write header
		fw.write("@RELATION " + FilenameUtils.getBaseName(arffFile.getAbsolutePath()) + "\n\n");
		for (final DistanceMeasure distance : getDistances()) {
			fw.write("@ATTRIBUTE " + distance.name() + " REAL\n");
		}
		if (writeClass) {
			fw.write("@ATTRIBUTE class {");
			final StringBuilder sb2 = new StringBuilder();
			for (final ClassLabel classLabel : ClassLabel.valuesArray()) {
				if (!"".equals(sb2.toString())) {
					sb2.append(",");
				}
				sb2.append(classLabel.name());
			}
			fw.write(sb2.toString() + "}\n");
		}
		fw.write("\n");
		fw.write("@DATA\n");
		// write data

		for (final String proteinPair : proteinPairKeysList) {

			final StringBuilder sb = new StringBuilder();
			for (final DistanceMeasure distance : getDistances()) {
				final double num = getValue(distance, proteinPair);
				if (!"".equals(sb.toString())) {
					sb.append(",");
				}
				if (!Double.isNaN(num)) {
					sb.append(num);
				} else {
					sb.append("?");
				}
			}
			if (writeClass) {
				final ClassLabel classLabel = trueClassifier.getClassLabel(
						getIndividualProteinsFromProteinGroupKey(getFirstProtein(proteinPair)),
						getIndividualProteinsFromProteinGroupKey(getSecondProtein(proteinPair)));
				if (classLabel == null) {
					trueClassifier.getClassLabel(getIndividualProteinsFromProteinGroupKey(getFirstProtein(proteinPair)),
							getIndividualProteinsFromProteinGroupKey(getSecondProtein(proteinPair)));
				}
				sb.append(",");

				sb.append(classLabel.name());
			}
			fw.write(sb.toString() + "\n");
		}

		fw.close();
		log.info("ARFF file created at " + arffFile.getAbsolutePath());

	}

	public String getFirstProtein(String proteinPair) {
		return proteinPair.split(PAIR_SEPARATOR)[0];
	}

	public String getSecondProtein(String proteinPair) {
		return proteinPair.split(PAIR_SEPARATOR)[1];
	}

	public List<String> getProteinPairs() {
		return proteinPairKeysList;
	}

	public ClassLabel getTrueClass(String proteinPair) throws IOException {
		final ClassLabel classLabel = trueClassifier.getClassLabel(
				getIndividualProteinsFromProteinGroupKey(getFirstProtein(proteinPair)),
				getIndividualProteinsFromProteinGroupKey(getSecondProtein(proteinPair)));
		return classLabel;
	}

	public TrueClassifier getTrueClassifier() {
		return trueClassifier;
	}
}
