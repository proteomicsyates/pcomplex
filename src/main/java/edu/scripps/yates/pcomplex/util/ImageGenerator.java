package edu.scripps.yates.pcomplex.util;

import java.awt.BasicStroke;
import java.io.File;
import java.io.IOException;

import org.apache.log4j.Logger;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;

public class ImageGenerator {
	private final static Logger log = Logger.getLogger(ImageGenerator.class);
	private static int width = 700;
	private static int height = 700;

	public static void generateImage(File imageFile, double[][] data, String title, String xAxis, String yAxis)
			throws IOException {

		final NumberAxis numberaxis_x = new NumberAxis(xAxis);
		numberaxis_x.setRange(0.0, 1.0);
		numberaxis_x.setTickUnit(new NumberTickUnit(0.1));
		final NumberAxis numberaxis_y = new NumberAxis(yAxis);
		numberaxis_y.setRange(0.0, 1.0);
		numberaxis_y.setTickUnit(new NumberTickUnit(0.1));
		final XYPlot xyplot = new XYPlot();
		xyplot.setDomainAxis(numberaxis_x);
		xyplot.setRangeAxis(numberaxis_y);

		// final XYSplineRenderer renderer = new XYSplineRenderer();
		final StandardXYItemRenderer renderer = new StandardXYItemRenderer();
		renderer.setSeriesStroke(1, new BasicStroke(0.0f));
		// renderer.setOutlineStrokeTableActive(false);
		xyplot.setRenderer(renderer);

		xyplot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
		final int numDataSet = 1;

		// experimental data
		final IntervalXYDataset intervalxydataset = getSeries(title, data);
		xyplot.setDataset(intervalxydataset);
		// renderer.setSeriesShapesVisible(0, false);
		// renderer.setSeriesLinesVisible(0, true);
		final JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, xyplot, true);
		final LegendTitle legend = chart.getLegend();
		legend.setPosition(RectangleEdge.BOTTOM);

		ChartUtils.saveChartAsPNG(imageFile, chart, width, height);
		log.info("Image saved at: " + imageFile.getAbsolutePath());
	}

	private static XYSeriesCollection getSeries(String title, double[][] data) {
		final ProgressCounter counter = new ProgressCounter(data.length, ProgressPrintingType.PERCENTAGE_STEPS, 0);
		counter.setShowRemainingTime(true);
		counter.setSuffix("creating image");
		final XYSeries xyseries = new XYSeries(title);
		for (int i = 0; i < data.length; i++) {
			counter.increment();
			final String printIfNecessary = counter.printIfNecessary();
			if (!"".equals(printIfNecessary)) {
				log.info(printIfNecessary);
			}
			final double x = data[i][0];
			final double y = data[i][1];
			xyseries.add(x, y);
		}
		return new XYSeriesCollection(xyseries);
	}
}
