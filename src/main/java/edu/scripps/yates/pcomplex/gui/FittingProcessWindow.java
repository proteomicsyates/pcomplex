package edu.scripps.yates.pcomplex.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.text.DecimalFormat;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import edu.scripps.yates.pcomplex.cofractionation.fitting.ModelPlot;
import edu.scripps.yates.utilities.dates.DatesUtil;

public class FittingProcessWindow extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5646969589172487789L;
	private final JPanel contentPane;
	private final JPanel northPanel;
	private final JLabel proteinACC;
	private final JLabel numGaussiansLabel;
	private final JLabel fitnessValueLabel;
	private final JPanel graphContainerPanel;
	private final JLabel previousFitnessValueLabel;
	private final DecimalFormat formatter = new DecimalFormat("#.##");
	private final JLabel lblR;
	private final JLabel rSquareLabel;
	private final JLabel lblAvgTimeFitting;
	private final JLabel avgTimeFittingLabel;

	/**
	 * Create the frame.
	 */
	public FittingProcessWindow() {
		setTitle("Chromatogram Profile");
		setBounds(100, 100, 450, 300);
		contentPane = new JPanel();
		contentPane.setPreferredSize(new Dimension(1600, 800));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);

		northPanel = new JPanel();
		northPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		northPanel.setBorder(new TitledBorder(null, "Status", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		contentPane.add(northPanel, BorderLayout.NORTH);
		final GridBagLayout gbl_northPanel = new GridBagLayout();
		gbl_northPanel.columnWidths = new int[] { 50, 50, 0, 0 };
		gbl_northPanel.rowHeights = new int[] { 0, 0, 0, 0 };
		gbl_northPanel.columnWeights = new double[] { 0.0, 0.0, 0.0, 0.0 };
		gbl_northPanel.rowWeights = new double[] { 0.0, 0.0, 1.0, Double.MIN_VALUE };
		northPanel.setLayout(gbl_northPanel);

		final JLabel lblProteins = new JLabel("Protein(s):");
		final GridBagConstraints gbc_lblProteins = new GridBagConstraints();
		gbc_lblProteins.anchor = GridBagConstraints.EAST;
		gbc_lblProteins.insets = new Insets(0, 0, 5, 5);
		gbc_lblProteins.gridx = 0;
		gbc_lblProteins.gridy = 0;
		northPanel.add(lblProteins, gbc_lblProteins);

		proteinACC = new JLabel("");
		final GridBagConstraints gbc_proteinACC = new GridBagConstraints();
		gbc_proteinACC.insets = new Insets(0, 0, 5, 5);
		gbc_proteinACC.anchor = GridBagConstraints.WEST;
		gbc_proteinACC.gridx = 1;
		gbc_proteinACC.gridy = 0;
		northPanel.add(proteinACC, gbc_proteinACC);

		final JLabel lblNewLabel_1 = new JLabel("Fitness value:");
		final GridBagConstraints gbc_lblNewLabel_1 = new GridBagConstraints();
		gbc_lblNewLabel_1.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel_1.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel_1.gridx = 2;
		gbc_lblNewLabel_1.gridy = 0;
		northPanel.add(lblNewLabel_1, gbc_lblNewLabel_1);

		final JLabel lblNewLabel = new JLabel("Number of gaussians used:");
		final GridBagConstraints gbc_lblNewLabel = new GridBagConstraints();
		gbc_lblNewLabel.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel.gridx = 0;
		gbc_lblNewLabel.gridy = 1;
		northPanel.add(lblNewLabel, gbc_lblNewLabel);

		numGaussiansLabel = new JLabel("");
		final GridBagConstraints gbc_numGaussiansLabel = new GridBagConstraints();
		gbc_numGaussiansLabel.insets = new Insets(0, 0, 5, 5);
		gbc_numGaussiansLabel.gridx = 1;
		gbc_numGaussiansLabel.gridy = 1;
		northPanel.add(numGaussiansLabel, gbc_numGaussiansLabel);

		fitnessValueLabel = new JLabel("");
		final GridBagConstraints gbc_errorLabel = new GridBagConstraints();
		gbc_errorLabel.insets = new Insets(0, 0, 5, 0);
		gbc_errorLabel.gridx = 3;
		gbc_errorLabel.gridy = 0;
		northPanel.add(fitnessValueLabel, gbc_errorLabel);

		final JLabel lblPreviousR = new JLabel("Previous fitness value:");
		final GridBagConstraints gbc_lblPreviousR = new GridBagConstraints();
		gbc_lblPreviousR.anchor = GridBagConstraints.EAST;
		gbc_lblPreviousR.insets = new Insets(0, 0, 5, 5);
		gbc_lblPreviousR.gridx = 2;
		gbc_lblPreviousR.gridy = 1;
		northPanel.add(lblPreviousR, gbc_lblPreviousR);

		previousFitnessValueLabel = new JLabel("");
		final GridBagConstraints gbc_lblPreviousrsquaredlabel = new GridBagConstraints();
		gbc_lblPreviousrsquaredlabel.insets = new Insets(0, 0, 5, 0);
		gbc_lblPreviousrsquaredlabel.gridx = 3;
		gbc_lblPreviousrsquaredlabel.gridy = 1;
		northPanel.add(previousFitnessValueLabel, gbc_lblPreviousrsquaredlabel);

		lblAvgTimeFitting = new JLabel("Avg time fitting:");
		final GridBagConstraints gbc_lblAvgTimeFitting = new GridBagConstraints();
		gbc_lblAvgTimeFitting.anchor = GridBagConstraints.EAST;
		gbc_lblAvgTimeFitting.insets = new Insets(0, 0, 0, 5);
		gbc_lblAvgTimeFitting.gridx = 0;
		gbc_lblAvgTimeFitting.gridy = 2;
		northPanel.add(lblAvgTimeFitting, gbc_lblAvgTimeFitting);

		avgTimeFittingLabel = new JLabel("");
		final GridBagConstraints gbc_avgTimeFittingLabel = new GridBagConstraints();
		gbc_avgTimeFittingLabel.anchor = GridBagConstraints.NORTH;
		gbc_avgTimeFittingLabel.insets = new Insets(0, 0, 0, 5);
		gbc_avgTimeFittingLabel.gridx = 1;
		gbc_avgTimeFittingLabel.gridy = 2;
		northPanel.add(avgTimeFittingLabel, gbc_avgTimeFittingLabel);

		lblR = new JLabel("R^2");
		final GridBagConstraints gbc_lblR = new GridBagConstraints();
		gbc_lblR.anchor = GridBagConstraints.EAST;
		gbc_lblR.insets = new Insets(0, 0, 0, 5);
		gbc_lblR.gridx = 2;
		gbc_lblR.gridy = 2;
		northPanel.add(lblR, gbc_lblR);

		rSquareLabel = new JLabel("");
		final GridBagConstraints gbc_rSquareLabel = new GridBagConstraints();
		gbc_rSquareLabel.anchor = GridBagConstraints.EAST;
		gbc_rSquareLabel.gridx = 3;
		gbc_rSquareLabel.gridy = 2;
		northPanel.add(rSquareLabel, gbc_rSquareLabel);

		graphContainerPanel = new JPanel();
		contentPane.add(graphContainerPanel, BorderLayout.CENTER);
		graphContainerPanel.setLayout(new GridLayout(4, 1, 0, 10));

		pack();
		final java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
		final java.awt.Dimension dialogSize = getSize();
		setLocation((screenSize.width - dialogSize.width) / 2, (screenSize.height - dialogSize.height) / 2);
		setVisible(true);
	}

	public void showGraph(ModelPlot modelPlot, Double previousRSquared) {
		fitnessValueLabel.setText(formatter.format(modelPlot.getFitnessValue()));
		proteinACC.setText(modelPlot.getAcc());
		numGaussiansLabel.setText(String.valueOf(modelPlot.getFitModel().getFittedGaussians().size()));
		if (previousRSquared != null) {
			previousFitnessValueLabel.setText(formatter.format(previousRSquared));
		} else {
			previousFitnessValueLabel.setText("-");
		}
		rSquareLabel.setText(formatter.format(modelPlot.getrSquare()));
		graphContainerPanel.removeAll();
		// spc
		final JFreeChart chart4 = modelPlot.getChart("Spectral counts profile", true, false, false, false, false);
		final ChartPanel graphPanel4 = new ChartPanel(chart4);
		graphPanel4.setFillZoomRectangle(false);
		graphPanel4.setBorder(new TitledBorder(null, "Chart", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		graphContainerPanel.add(graphPanel4);

		// experimental raw profile
		// final JFreeChart chart0 = modelPlot.getChart("Experimental RAW
		// profile", false, true, false, false, false);
		// final ChartPanel graphPanel0 = new ChartPanel(chart0);
		// graphPanel0.setFillZoomRectangle(false);
		// graphPanel0.setBorder(new TitledBorder(null, "Chart",
		// TitledBorder.LEADING, TitledBorder.TOP, null, null));
		// graphContainerPanel.add(graphPanel0);

		// experimental profile with fitted gaussians
		final JFreeChart chart = modelPlot.getChart("Experimental processed data", false, false, true, true, false);
		final ChartPanel graphPanel = new ChartPanel(chart);
		graphPanel.setFillZoomRectangle(false);
		graphPanel.setBorder(new TitledBorder(null, "Chart", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		graphContainerPanel.add(graphPanel);

		// // fitted gaussians
		// final JFreeChart chart2 = modelPlot.getChart("Fitted gaussians",
		// false, false, false, true, false);
		// final ChartPanel graphPanel2 = new ChartPanel(chart2);
		// graphPanel2.setFillZoomRectangle(false);
		// graphPanel2.setBorder(new TitledBorder(null, "Chart",
		// TitledBorder.LEADING, TitledBorder.TOP, null, null));
		// graphContainerPanel.add(graphPanel2);

		// total fitted gaussians
		final JFreeChart chart3 = modelPlot.getChart("Total fitted gaussians", false, false, true, false, true);
		final ChartPanel graphPanel3 = new ChartPanel(chart3);
		graphPanel3.setFillZoomRectangle(false);
		graphPanel3.setBorder(new TitledBorder(null, "Chart", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		graphContainerPanel.add(graphPanel3);

		setVisible(true);
	}

	public void setAverageFittingTime(double mean) {
		final String timeString = DatesUtil.getDescriptiveTimeFromMillisecs(mean);
		avgTimeFittingLabel.setText(timeString);
	}
}
