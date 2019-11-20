package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import gnu.trove.list.TDoubleList;
import gnu.trove.map.hash.THashMap;

public class FittingCache {
	private final static Map<SeparationExperiment, FittingCache> map = new THashMap<SeparationExperiment, FittingCache>();
	private static final String PROTEIN = "P";
	private static final String GAUSSIAN = "G";
	private final SeparationExperiment exp;
	private final Map<String, List<MyGaussianFit>> fitsByProtein = new THashMap<String, List<MyGaussianFit>>();

	private FittingCache(SeparationExperiment exp) {
		this.exp = exp;
		loadFromFile();
	}

	private void loadFromFile() {
		final File file = getFile();
		if (!file.exists() || file.length() == 0l) {
			return;
		}
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(file));
			String line = reader.readLine().trim();
			String acc = null;
			while (line != null) {
				try {
					line = line.trim();
					if (!"".equals(line)) {
						final String[] split = line.split("\t");
						if (line.startsWith(PROTEIN)) {
							acc = split[1].trim();
							final List<MyGaussianFit> fittings = new ArrayList<MyGaussianFit>();
							fitsByProtein.put(acc, fittings);
						} else if (line.startsWith(GAUSSIAN)) {
							final double h = Double.valueOf(split[1]);
							final double mean = Double.valueOf(split[2]);
							final double sigma = Double.valueOf(split[3]);
							final MyGaussianFit fit = new MyGaussianFit(h, mean, sigma);
							fitsByProtein.get(acc).add(fit);
						}
					}
				} finally {
					line = reader.readLine();
				}
			}
			reader.close();
		} catch (final IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException("Error loading file " + getFile().getAbsolutePath());
		}
	}

	private File getFile() {
		return new File(exp.getFile().getParentFile() + File.separator + exp.getProjectName() + "_apex_fittings.txt");
	}

	public static FittingCache getInstance(SeparationExperiment exp) {
		if (!map.containsKey(exp)) {
			map.put(exp, new FittingCache(exp));
		}
		return map.get(exp);
	}

	public List<MyGaussianFit> getFittingsForProtein(String acc) {
		return fitsByProtein.get(acc);
	}

	public boolean containsFittingsForProtein(String acc) {
		return fitsByProtein.containsKey(acc);
	}

	public void put(String acc, List<MyGaussianFit> fits) throws IOException {

		writeToFile(acc, fits);

		fitsByProtein.put(acc, fits);
	}

	private void writeToFile(String acc, List<MyGaussianFit> fits) throws IOException {
		FileWriter fw = null;
		try {
			fw = new FileWriter(getFile(), true);

			fw.write(PROTEIN + "\t" + acc + "\n");

			for (final MyGaussianFit myGaussianFit : fits) {
				final TDoubleList fittedParameters = myGaussianFit.getFittedParameters();
				fw.write(GAUSSIAN + "\t" + FittingUtil.getHeight(fittedParameters) + "\t"
						+ FittingUtil.getMean(fittedParameters) + "\t" + FittingUtil.getSigma(fittedParameters) + "\n");
			}

		} finally {
			if (fw != null) {
				fw.close();
			}
		}

	}

}
