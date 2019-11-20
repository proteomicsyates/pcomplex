package edu.scripps.yates.pcomplex.cofractionation.fitting;

import java.awt.Color;
import java.util.List;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jgap.Population;

import edu.scripps.yates.pcomplex.cofractionation.MyGaussianFit;
import gnu.trove.list.TDoubleList;

/**
 * class overlay the fitted model over the distribution
 *
 * @author diego and Salva
 * 
 */
public class ModelPlot {
	private final List<double[]> profile;
	private final AbstractMultipleFitModel fitModel;
	private boolean _final;
	private final String acc;
	private final List<double[]> rawProfile;
	private final List<int[]> spcProfile;
	private final Population population;

	/**
	 * @param population
	 * @param list
	 * @param title
	 */
	public ModelPlot(String acc, List<int[]> spcProfile, List<double[]> rawProfile, List<double[]> histdata,
			AbstractMultipleFitModel fitModel, Population population) {
		profile = histdata;
		this.spcProfile = spcProfile;
		this.rawProfile = rawProfile;
		this.fitModel = fitModel;
		this.acc = acc;
		this.population = population;
	}

	public JFreeChart getChart(String title, boolean showSPC, boolean showRawExperimentalData,
			boolean showExperimentalProcessedData, boolean showFittedGaussians, boolean showTotal) {

		final NumberAxis numberaxis_x = new NumberAxis("Fraction");
		String yLabel = "Norm SPC";
		if (showSPC) {
			yLabel = "SPC";
		}
		final NumberAxis numberaxis_y = new NumberAxis(yLabel);

		final XYPlot xyplot = new XYPlot();
		xyplot.setDomainAxis(numberaxis_x);
		xyplot.setRangeAxis(numberaxis_y);
		if (false) {
			// this is to see vertical lines in each point
			final XYBarRenderer xybarrenderer = new XYBarRenderer();
			xybarrenderer.setShadowVisible(false);
			xybarrenderer.setMargin(0.95);
			xyplot.setRenderer(xybarrenderer);
		} else {
			xyplot.setRenderer(new XYLineAndShapeRenderer());
		}
		xyplot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

		int numDataSet = 0;
		if (showSPC) {
			// experimental data
			final IntervalXYDataset intervalxydataset = createDatasetFromSPCProfile(spcProfile);
			xyplot.setDataset(numDataSet++, intervalxydataset);
			final XYBarRenderer xybarrenderer = new XYBarRenderer();
			xybarrenderer.setShadowVisible(false);
			xybarrenderer.setMargin(0.95);
			xyplot.setRenderer(numDataSet - 1, xybarrenderer);
		}
		if (showRawExperimentalData) {
			// experimental data
			final IntervalXYDataset intervalxydataset = createDatasetFromProfile(rawProfile);

			xyplot.setDataset(numDataSet++, intervalxydataset);
		}

		if (showExperimentalProcessedData) {
			// experimental data
			final IntervalXYDataset intervalxydataset = createDatasetFromProfile(profile);

			xyplot.setDataset(numDataSet++, intervalxydataset);
		}
		// renderer.setSeriesPaint(0, Color.BLUE);
		// renderer.setSeriesPaint(1, Color.RED);
		// renderer.setSeriesPaint(2, Color.BLACK);

		if (showFittedGaussians || showTotal) {
			if (fitModel.getFittedGaussians() != null && !fitModel.getFittedGaussians().isEmpty()) {
				// individual gaussians
				if (showFittedGaussians) {
					for (final MyGaussianFit gaussian : fitModel.getFittedGaussians()) {

						final XYDataset xydataset2 = createXYDatasetFromNormalDistribution(gaussian,
								"Gaussian " + numDataSet);
						xyplot.setDataset(numDataSet, xydataset2);
						xyplot.setRenderer(numDataSet, new StandardXYItemRenderer());

						numDataSet++;
					}
				}
				if (showTotal) {
					// total
					final XYDataset xydataset3 = createXYDatasetFromNormalDistribution(fitModel.getFittedGaussians(),
							"TOTAL");
					xyplot.setDataset(numDataSet, xydataset3);
					final StandardXYItemRenderer renderer = new StandardXYItemRenderer();
					renderer.setSeriesPaint(numDataSet, Color.BLACK);
					xyplot.setRenderer(numDataSet, renderer);
				}
			}
		}
		return new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, xyplot, true);
	}

	private XYDataset createXYDatasetFromNormalDistribution(MyGaussianFit gaussian, String seriesName) {

		final XYSeries xyseries = new XYSeries(seriesName);

		for (int j = 0; j < profile.size(); j++) {
			xyseries.add(profile.get(j)[0], gaussian.getGaussianY(profile.get(j)[0]));
		}

		return new XYSeriesCollection(xyseries);

	}

	private XYDataset createXYDatasetFromNormalDistribution(List<MyGaussianFit> gaussians, String seriesName) {

		final XYSeries xyseries = new XYSeries(seriesName);

		for (int j = 0; j < profile.size(); j++) {
			double yValue = 0.0;
			final double x = profile.get(j)[0];
			for (final MyGaussianFit gaussian : gaussians) {
				yValue += gaussian.getGaussianY(x);
			}
			xyseries.add(x, yValue);
		}

		return new XYSeriesCollection(xyseries);

	}

	private IntervalXYDataset createDatasetFromProfile(List<double[]> profile) {
		// log.info("Histogram values are:");
		// StringBuilder sb = new StringBuilder();
		// for (double[] values : experimentalHistogramData) {
		// sb.append("(" + values[0] + "," + values[1] + "), ");
		// }
		// // log.info(sb.toString());

		final XYSeries xyseries = new XYSeries("Experimental data");
		for (int i = 0; i < profile.size(); i++) {
			final double x = profile.get(i)[0];
			final double y = profile.get(i)[1];
			xyseries.add(x, y);
		}
		return new XYSeriesCollection(xyseries);
	}

	private IntervalXYDataset createDatasetFromSPCProfile(List<int[]> profile) {
		// log.info("Histogram values are:");
		// StringBuilder sb = new StringBuilder();
		// for (double[] values : experimentalHistogramData) {
		// sb.append("(" + values[0] + "," + values[1] + "), ");
		// }
		// // log.info(sb.toString());

		final XYSeries xyseries = new XYSeries("Spectral counts");
		for (int i = 0; i < profile.size(); i++) {
			final int x = profile.get(i)[0];
			final int y = profile.get(i)[1];
			xyseries.add(x, y);
		}
		return new XYSeriesCollection(xyseries);
	}

	/**
	 * @return the fitModelSimple
	 */
	public AbstractMultipleFitModel getFitModel() {
		return fitModel;
	}

	public int getNumBins() {
		return profile.size();
	}

	public double getrSquare() {
		return fitModel.getRSquare();
	}

	public boolean isFinal() {
		return _final;
	}

	public void setFinal(boolean _final) {
		this._final = _final;
	}

	public TDoubleList getProcessedProfile() {
		return fitModel.getProcessedProfile();
	}

	public String getAcc() {
		return acc;
	}

	public double getFitnessValue() {
		return fitModel.getFitnessValue();
	}

	public Population getPopulation() {
		return population;
	}

}
