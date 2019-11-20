package edu.scripps.yates.pcomplex.epic;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

/**
 * It reads a Summary file that was copied and pasted in a txt file from the web
 * browser output of EPIC:<br>
 * Average size of predicted complexes is: 9.64814814815 <br>
 * mmr 0.217029387745 <br>
 * overlapp 0.311111 <br>
 * simcoe 0.366667 <br>
 * mean_simcoe_overlap 0.359259 <br>
 * sensetivity 0.730861819933 <br>
 * ppv 0.229549957543 <br>
 * accuracy 0.409596508451 <br>
 * sep 0.161019014645 <br>
 * composite score is: 0.937736896196
 * 
 * @author salvador
 *
 */
public class EpicSummary {
	private double averageSize;
	private double mmr;
	private double overlap;
	private double simcoe;
	private double meanSimcoeOverlap;
	private double sensitivity;
	private double ppv;
	private double accuracy;
	private double sep;
	private double compositeScore;

	public static String getHeaderLine() {
		return "averageSize\tmmr\toverlap\tsimcoe\tmeanSimcoeOverlap\tsensitivity\tppv\taccuracy\tsep\tcompositeScore";
	}

	public String getString() {
		final StringBuilder sb = new StringBuilder();
		sb.append(averageSize).append("\t").append(mmr).append("\t").append(overlap).append("\t").append(simcoe)
				.append("\t").append(meanSimcoeOverlap).append("\t").append(sensitivity).append("\t").append(ppv)
				.append("\t").append(accuracy).append("\t").append(sep).append("\t").append(compositeScore);
		return sb.toString();
	}

	public EpicSummary(File epicSummaryFile) throws IOException {
		final List<String> lines = Files.readAllLines(epicSummaryFile.toPath());
		for (final String line : lines) {
			if (!line.contains("is:")) {
				final String[] split = line.split("\t");
				final double number = Double.valueOf(split[1]);
				if ("mmr".equals(split[0])) {
					mmr = number;
				} else if ("overlapp".equals(split[0])) {
					overlap = number;
				} else if ("simcoe".equals(split[0])) {
					simcoe = number;
				} else if ("mean_simcoe_overlap".equals(split[0])) {
					meanSimcoeOverlap = number;
				} else if ("sensetivity".equals(split[0])) {
					sensitivity = number;
				} else if ("ppv".equals(split[0])) {
					ppv = number;
				} else if ("accuracy".equals(split[0])) {
					accuracy = number;
				} else if ("sep".equals(split[0])) {
					sep = number;
				}
			} else {
				final String[] split = line.split("\\s");
				final Double number = Double.valueOf(split[split.length - 1]);
				if ("composite".equalsIgnoreCase(split[0])) {
					compositeScore = number;
				} else if ("average".equalsIgnoreCase(split[0])) {
					averageSize = number;
				}
			}
		}
	}

	public double getAverageSize() {
		return averageSize;
	}

	public double getMmr() {
		return mmr;
	}

	public double getOverlap() {
		return overlap;
	}

	public double getSimcoe() {
		return simcoe;
	}

	public double getMeanSimcoeOverlap() {
		return meanSimcoeOverlap;
	}

	public double getSensitivity() {
		return sensitivity;
	}

	public double getPpv() {
		return ppv;
	}

	public double getAccuracy() {
		return accuracy;
	}

	public double getSep() {
		return sep;
	}

	public double getCompositeScore() {
		return compositeScore;
	}

}
