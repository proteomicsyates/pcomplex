package edu.scripps.yates.pcomplex.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import edu.scripps.yates.pcomplex.ProteinComplexAnalyzer;
import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.epic.EpicResultComparator;
import edu.scripps.yates.pcomplex.model.Fraction;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;
import edu.scripps.yates.pcomplex.model.SeparationExperiment;
import edu.scripps.yates.pcomplex.util.DataType;
import edu.scripps.yates.pcomplex.util.PComplexUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotGeneMapping;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class ComplexesWindow extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5646969589172487789L;
	private final JPanel contentPane;
	private final JPanel graphContainerPanel;
	private static final String ACCS = "accs";
	private static final String GENES = "genes";
	private final JPanel panel;
	private final JButton btnSubmitProfileData;
	private final JButton btnSubmitProteinComplexes;
	private final JPanel panel_1;
	private final JPanel panel_2;
	private File currentDirectoryProfiles = new File(
			"Z:\\share\\Salva\\data\\asaviola\\protein complexes\\experiments");
	private File currentDirectoryComplexes = new File(
			"Z:\\share\\Salva\\data\\asaviola\\protein complexes\\ML\\matrixes");
	private SeparationExperiment experiment;
	private final DefaultListModel<ProteinComplex> defaultProteinComplexesModel;
	private static final Logger log = Logger.getLogger(ComplexesWindow.class);
	private final JComboBox<DataType> comboBoxDataType;
	private final JLabel lblDataType;
	private final JCheckBox chckbxSeparated;
	private final JList<ProteinComplex> complexList;
	private final JPanel panel_3;
	private final JPanel panel_4;
	private final JCheckBox chckbxEnable;
	private final JPanel panel_5;
	private final JLabel lblSmoothWindow;
	private final JTextField smoothWindowText;
	private final JPanel panel_6;
	private final JTextField searchProteinText;
	private final JLabel numberComplexesLabel;
	private final JPanel panel_7;
	private final JLabel proteinsInSelectedClusterLabel;
	private final JPanel panel_8;
	private final JLabel lblNewLabel;
	private final JComboBox<String> accsOrGenesComboBox;
	private final JPanel panel_9;
	private final JPanel verticalRightPanel;
	private final JPanel panel_11;
	private final JLabel lblNewLabel_1;
	private final JButton btnNewButton;
	private final List<ReferenceDBPanel> referenceDBPanels = new ArrayList<ReferenceDBPanel>();

	public static void main(String[] args) {
		ProteinComplexAnalyzer.useComplexPortalDB = true;
		ProteinComplexAnalyzer.useCoreCorumDB = true;
		ProteinComplexAnalyzer.useHUMAP = true;
		final ComplexesWindow window = new ComplexesWindow();
		window.setVisible(true);
	}

	/**
	 * Create the frame.
	 */
	public ComplexesWindow() {
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (UnsupportedLookAndFeelException | ClassNotFoundException | InstantiationException
				| IllegalAccessException e) {
			e.printStackTrace();
		}
		setTitle("Chromatogram Profile");
		setBounds(100, 100, 450, 300);
		contentPane = new JPanel();
		contentPane.setPreferredSize(new Dimension(1600, 800));
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		contentPane.setLayout(new BorderLayout(0, 0));
		setContentPane(contentPane);

		panel = new JPanel();
		panel.setPreferredSize(new Dimension(250, 800));
		panel.setSize(new Dimension(250, 0));
		panel.setMinimumSize(new Dimension(250, 10));
		contentPane.add(panel, BorderLayout.WEST);
		panel.setLayout(new BorderLayout(0, 20));

		panel_1 = new JPanel();
		panel.add(panel_1, BorderLayout.NORTH);
		final GridBagLayout gbl_panel_1 = new GridBagLayout();
		gbl_panel_1.columnWidths = new int[] { 221, 0 };
		gbl_panel_1.rowHeights = new int[] { 0, 0, 0, 23, 0, 0 };
		gbl_panel_1.columnWeights = new double[] { 1.0, Double.MIN_VALUE };
		gbl_panel_1.rowWeights = new double[] { 0.0, 1.0, 0.0, 1.0, 1.0, Double.MIN_VALUE };
		panel_1.setLayout(gbl_panel_1);

		lblDataType = new JLabel("Data type:");
		lblDataType.setHorizontalAlignment(SwingConstants.RIGHT);

		final JPanel panelHorizontal = new JPanel();
		panelHorizontal.setLayout(new GridLayout(1, 2, 0, 0));
		panelHorizontal.add(lblDataType);

		comboBoxDataType = new JComboBox<DataType>();
		comboBoxDataType.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				dataTypeSelected();
			}
		});
		final GridBagConstraints gbc_comboBoxDataType = new GridBagConstraints();
		gbc_comboBoxDataType.insets = new Insets(0, 0, 5, 0);
		gbc_comboBoxDataType.fill = GridBagConstraints.HORIZONTAL;
		gbc_comboBoxDataType.gridx = 0;
		gbc_comboBoxDataType.gridy = 0;
		panelHorizontal.add(comboBoxDataType);
		panel_1.add(panelHorizontal, gbc_comboBoxDataType);

		chckbxSeparated = new JCheckBox("Separate graphs");
		chckbxSeparated.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				separatedSelected();
			}
		});

		panel_8 = new JPanel();
		final GridBagConstraints gbc_panel_8 = new GridBagConstraints();
		gbc_panel_8.insets = new Insets(0, 0, 5, 0);
		gbc_panel_8.fill = GridBagConstraints.BOTH;
		gbc_panel_8.gridx = 0;
		gbc_panel_8.gridy = 1;
		panel_1.add(panel_8, gbc_panel_8);
		panel_8.setLayout(new GridLayout(1, 2, 0, 0));

		lblNewLabel = new JLabel("ACCs or Genes:");
		lblNewLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		panel_8.add(lblNewLabel);

		accsOrGenesComboBox = new JComboBox<String>();
		accsOrGenesComboBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				accsorGenesSelected();
			}
		});

		accsOrGenesComboBox.setModel(new DefaultComboBoxModel<String>(new String[] { ACCS, GENES }));
		panel_8.add(accsOrGenesComboBox);
		final GridBagConstraints gbc_chckbxNewCheckBox = new GridBagConstraints();
		gbc_chckbxNewCheckBox.insets = new Insets(0, 0, 5, 0);
		gbc_chckbxNewCheckBox.gridx = 0;
		gbc_chckbxNewCheckBox.gridy = 2;
		panel_1.add(chckbxSeparated, gbc_chckbxNewCheckBox);

		panel_4 = new JPanel();
		panel_4.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Elution profile processing",
				TitledBorder.LEADING, TitledBorder.TOP, null, new Color(0, 0, 0)));
		final GridBagConstraints gbc_panel_4 = new GridBagConstraints();
		gbc_panel_4.insets = new Insets(0, 0, 5, 0);
		gbc_panel_4.fill = GridBagConstraints.BOTH;
		gbc_panel_4.gridx = 0;
		gbc_panel_4.gridy = 3;
		panel_1.add(panel_4, gbc_panel_4);
		panel_4.setLayout(new GridLayout(2, 1, 0, 0));

		chckbxEnable = new JCheckBox("enable");
		chckbxEnable.setToolTipText("Check this to enable smooth of the elution profile");
		chckbxEnable.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final boolean enabled = chckbxEnable.isSelected();
				lblSmoothWindow.setEnabled(enabled);
				smoothWindowText.setEnabled(enabled);
				showSelectedComplex();
			}
		});
		panel_4.add(chckbxEnable);

		panel_5 = new JPanel();
		panel_4.add(panel_5);

		lblSmoothWindow = new JLabel("Smooth window:");
		lblSmoothWindow.setEnabled(false);
		panel_5.add(lblSmoothWindow);

		smoothWindowText = new JTextField();
		smoothWindowText.setToolTipText(
				"Number of fractions that would be used to smooth the values on each fraction by averaging them");
		smoothWindowText.setText("3");
		smoothWindowText.setEnabled(false);
		panel_5.add(smoothWindowText);
		smoothWindowText.setColumns(10);
		smoothWindowText.getDocument().addDocumentListener(new DocumentListener() {

			@Override
			public void removeUpdate(DocumentEvent e) {
				if (!"".equals(smoothWindowText.getText())) {
					showSelectedComplex();
				}

			}

			@Override
			public void insertUpdate(DocumentEvent e) {
				if (!"".equals(smoothWindowText.getText())) {
					showSelectedComplex();
				}
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				if (!"".equals(smoothWindowText.getText())) {
					showSelectedComplex();
				}
			}
		});
		panel_3 = new JPanel();
		panel_3.setBorder(new TitledBorder(null, "Input data", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		final GridBagConstraints gbc_panel_3 = new GridBagConstraints();
		gbc_panel_3.fill = GridBagConstraints.BOTH;
		gbc_panel_3.gridx = 0;
		gbc_panel_3.gridy = 4;
		panel_1.add(panel_3, gbc_panel_3);
		panel_3.setLayout(new GridLayout(2, 0, 0, 10));

		btnSubmitProfileData = new JButton("Submit profile data");
		panel_3.add(btnSubmitProfileData);

		btnSubmitProteinComplexes = new JButton("Submit protein complexes");
		panel_3.add(btnSubmitProteinComplexes);
		btnSubmitProteinComplexes.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				submitProteinComplexes();
			}
		});
		btnSubmitProfileData.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				addProfileData();
			}
		});

		panel_2 = new JPanel();
		panel.add(panel_2, BorderLayout.CENTER);
		panel_2.setLayout(new BorderLayout(0, 10));

		defaultProteinComplexesModel = new DefaultListModel<ProteinComplex>();
		final DefaultListModel<ProteinComplex> complexesModel = new DefaultListModel<ProteinComplex>();
		complexList = new JList<ProteinComplex>(complexesModel);
		complexList.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		complexList.setVisibleRowCount(50);
		complexList.addListSelectionListener(new ListSelectionListener() {

			@Override
			public void valueChanged(ListSelectionEvent e) {
				if (!e.getValueIsAdjusting()) {
					showSelectedComplex();

				}
			}
		});

		panel_6 = new JPanel();
		panel_2.add(panel_6, BorderLayout.NORTH);
		panel_6.setLayout(new GridLayout(2, 2, 5, 5));

		final JLabel lblSearch = new JLabel("Search:");
		lblSearch.setHorizontalAlignment(SwingConstants.RIGHT);
		panel_6.add(lblSearch);

		searchProteinText = new JTextField();
		searchProteinText.getDocument().addDocumentListener(new DocumentListener() {

			@Override
			public void removeUpdate(DocumentEvent e) {
				try {
					filterList(searchProteinText.getText());
				} catch (final IOException e1) {
					e1.printStackTrace();
				}
			}

			@Override
			public void insertUpdate(DocumentEvent e) {
				try {
					filterList(searchProteinText.getText());
				} catch (final IOException e1) {
					e1.printStackTrace();
				}
			}

			@Override
			public void changedUpdate(DocumentEvent e) {
				try {
					filterList(searchProteinText.getText());
				} catch (final IOException e1) {
					e1.printStackTrace();
				}
			}
		});
		panel_6.add(searchProteinText);
		searchProteinText.setColumns(10);

		final JLabel lblNumberOfClusters = new JLabel("Number of clusters:");
		lblNumberOfClusters.setHorizontalAlignment(SwingConstants.RIGHT);
		panel_6.add(lblNumberOfClusters);

		numberComplexesLabel = new JLabel("0");
		panel_6.add(numberComplexesLabel);

		panel_2.add(new JScrollPane(complexList));

		panel_7 = new JPanel();
		panel_2.add(panel_7, BorderLayout.SOUTH);
		panel_7.setLayout(new GridLayout(1, 2, 5, 0));

		final JLabel lblProteinsInSelected = new JLabel("Selected cluster:");
		lblProteinsInSelected.setHorizontalAlignment(SwingConstants.RIGHT);
		panel_7.add(lblProteinsInSelected);

		proteinsInSelectedClusterLabel = new JLabel("0 proteins");
		proteinsInSelectedClusterLabel.setHorizontalAlignment(SwingConstants.LEFT);
		panel_7.add(proteinsInSelectedClusterLabel);

		graphContainerPanel = new JPanel();
		contentPane.add(graphContainerPanel, BorderLayout.CENTER);
		final GridBagLayout gbl_graphContainerPanel = new GridBagLayout();
		gbl_graphContainerPanel.columnWidths = new int[] { 0, 0 };
		gbl_graphContainerPanel.rowHeights = new int[] { 0, 0 };
		gbl_graphContainerPanel.columnWeights = new double[] { 1.0, Double.MIN_VALUE };
		gbl_graphContainerPanel.rowWeights = new double[] { 1.0, Double.MIN_VALUE };
		graphContainerPanel.setLayout(gbl_graphContainerPanel);

		panel_9 = new JPanel();
		contentPane.add(panel_9, BorderLayout.EAST);

		verticalRightPanel = new JPanel();
		verticalRightPanel.setBorder(
				new TitledBorder(null, "Reference DB(s)", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_9.add(verticalRightPanel);
		final GridBagLayout gbl_panel_10 = new GridBagLayout();
		gbl_panel_10.columnWeights = new double[] { 0.0, Double.MIN_VALUE };
		gbl_panel_10.rowWeights = new double[] { 0.0, 0.0 };
		verticalRightPanel.setLayout(gbl_panel_10);

		panel_11 = new JPanel();
		final GridBagConstraints gbc_panel_11 = new GridBagConstraints();
		gbc_panel_11.fill = GridBagConstraints.BOTH;
		gbc_panel_11.insets = new Insets(0, 0, 5, 0);
		gbc_panel_11.gridx = 0;
		gbc_panel_11.gridy = 0;
		verticalRightPanel.add(panel_11, gbc_panel_11);

		lblNewLabel_1 = new JLabel("Load Reference DB");
		panel_11.add(lblNewLabel_1);

		btnNewButton = new JButton("Load");
		btnNewButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				selectReferenceDBFile();
			}
		});
		panel_11.add(btnNewButton);

		loadReferenceDBs();

		pack();
		final java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
		final java.awt.Dimension dialogSize = getSize();
		setLocation((screenSize.width - dialogSize.width) / 2, (screenSize.height - dialogSize.height) / 2);

		// comboBoxDataType.addItem(DataType.NSAF);
		comboBoxDataType.addItem(DataType.SPC);
	}

	private void loadReferenceDBs() {
		final ReferenceDBPanel referenceDBPanel = new ReferenceDBPanel();
		addReferenceDBPanel(referenceDBPanel);
		final List<ProteinComplexDB> dBs = ProteinComplexAnalyzer.getDBs();
		for (final ProteinComplexDB proteinComplexDB : dBs) {
			referenceDBPanel.addDB(proteinComplexDB);
		}
	}

	private void addReferenceDBPanel(ReferenceDBPanel referenceDBPanel) {
		referenceDBPanels.add(referenceDBPanel);

		final GridBagConstraints gbc_referenceDBPanel = new GridBagConstraints();
		gbc_referenceDBPanel.fill = GridBagConstraints.BOTH;
		gbc_referenceDBPanel.gridx = 0;
		gbc_referenceDBPanel.gridy = referenceDBPanels.size();
		verticalRightPanel.add(referenceDBPanel, gbc_referenceDBPanel);
	}

	protected void selectReferenceDBFile() {
		try {
			final JFileChooser fileChooser = new JFileChooser(currentDirectoryComplexes);
			fileChooser.showOpenDialog(this);
			final File selectedFile = fileChooser.getSelectedFile();

			loadReferenceDB(FilenameUtils.getBaseName(selectedFile.getAbsolutePath()), selectedFile);
		} catch (final Exception e) {
			showError(e);
		}
	}

	private void loadReferenceDB(String referenceDBName, File referenceDB) {
		// TODO Auto-generated method stub

	}

	protected void accsorGenesSelected() {
		if (accsOrGenesComboBox.getSelectedItem().equals(GENES)) {
			ProteinComplex.useGeneToPrint = true;
		} else {
			ProteinComplex.useGeneToPrint = false;
		}
		showSelectedComplex();
	}

	protected void filterList(String text) throws IOException {
		final DefaultListModel<ProteinComplex> listModel = (DefaultListModel<ProteinComplex>) complexList.getModel();
		try {
			listModel.clear();
		} catch (final Exception e) {
		}
		if ("".equals(text)) {
			for (int i = 0; i < defaultProteinComplexesModel.size(); i++) {
				listModel.addElement(defaultProteinComplexesModel.get(i));
			}
		} else {
			final ProteinComponent protein = new ProteinComponent(text, null);
			listModel.clear();
			for (int i = 0; i < defaultProteinComplexesModel.size(); i++) {
				if (defaultProteinComplexesModel.get(i).getComponentSet().contains(protein)) {
					listModel.addElement(defaultProteinComplexesModel.get(i));
				}
			}
		}
		updateNumberOfClusters();
	}

	protected void separatedSelected() {
		showSelectedComplex();

	}

	protected void showSelectedComplex() {
		setVisible(true);
		complexList.updateUI();
		if (complexList.isSelectionEmpty()) {
			return;
		}
		final ProteinComplex selectedComplex = complexList.getSelectedValue();
		// showError("Selected complex: " + selectedComplex, "Complex
		// selected:");
		showProfilesForComplex(selectedComplex);

		for (final ReferenceDBPanel referenceDBPanel : referenceDBPanels) {
			referenceDBPanel.loadOverlapsWithComplex(selectedComplex);
		}

	}

	protected void dataTypeSelected() {
		showSelectedComplex();
	}

	protected void submitProteinComplexes() {
		try {
			final JFileChooser fileChooser = new JFileChooser(currentDirectoryComplexes);
			fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
			fileChooser.showOpenDialog(this);
			final File selectedFile = fileChooser.getSelectedFile();
			currentDirectoryComplexes = selectedFile.getParentFile();

			List<ProteinComplex> proteincomplexes = null;
			try {
				proteincomplexes = loadProteinComplexes(selectedFile);
			} catch (final IllegalArgumentException e) {
				proteincomplexes = loadProteinComplexesFromEPICFolder(selectedFile);
			}
			loadProteinComplexesOnList(proteincomplexes);
		} catch (final Exception e) {
			showError(e);
		}
	}

	private void loadProteinComplexesOnList(List<ProteinComplex> proteincomplexes) {
		log.info("Loading protein complexes " + proteincomplexes.size());
		defaultProteinComplexesModel.clear();
		final DefaultListModel<ProteinComplex> model = (DefaultListModel<ProteinComplex>) complexList.getModel();
		model.clear();
		for (final ProteinComplex proteinComplex : proteincomplexes) {
			defaultProteinComplexesModel.addElement(proteinComplex);
			model.addElement(proteinComplex);
		}
		updateNumberOfClusters();
	}

	private void updateNumberOfClusters() {
		numberComplexesLabel.setText(String.valueOf(complexList.getModel().getSize()));
	}

	protected void showProfilesForComplex(ProteinComplex selectedComplex) {
		showGraph(selectedComplex, getSelectedDataType(), isSeparated(), getSmoothWindow());
	}

	private Integer getSmoothWindow() {
		if (chckbxEnable.isSelected()) {
			try {
				final String text = smoothWindowText.getText();

				return Integer.valueOf(text);

			} catch (final NumberFormatException e) {
				showError("Smooth window only allow integer positive numbers", "Error in smooth window");
			}
		}
		return null;
	}

	private boolean isSeparated() {
		return chckbxSeparated.isSelected();
	}

	private DataType getSelectedDataType() {
		return (DataType) comboBoxDataType.getSelectedItem();
	}

	private List<ProteinComplex> loadProteinComplexes(File selectedFile) throws IOException {
		if (!selectedFile.isFile()) {
			throw new IllegalArgumentException("This is not a file!");
		}
		final List<String> readAllLines = Files.readAllLines(selectedFile.toPath());
		final List<ProteinComplex> ret = new ArrayList<ProteinComplex>();
		final TObjectIntHashMap<String> indexByHeader = new TObjectIntHashMap<String>();
		for (final String line : readAllLines) {
			if (line.startsWith("Type")) {
				int index = 0;
				final String[] split = line.split("\t");
				for (final String header : split) {
					indexByHeader.put(header, index);
					index++;
				}
				continue;
			}
			if (indexByHeader.isEmpty()) {
				throw new IllegalArgumentException("This is not the proper format. Try with reading EPIC results");
			}
			final ProteinComplex complex = new ProteinComplex(null);
			ret.add(complex);
			final List<String> components = new ArrayList<String>();
			if (!indexByHeader.isEmpty() && indexByHeader.containsKey("Components")) {
				if (line.contains("\t")) {
					final String[] split2 = line.split("\t");
					for (int i = indexByHeader.get("Components"); i < split2.length; i++) {
						components.add(split2[i]);
					}
				}
			} else {
				if (line.contains("\t")) {
					final String[] split = line.split("\t");
					for (final String acc : split) {
						components.add(acc);
					}
				} else {
					components.add(line);
				}
			}
			final String uniProtACC = FastaParser.getUniProtACC(components.get(0));
			if (uniProtACC != null) {
				for (final String acc : components) {
					String gene = null;
					final Map<String, Set<String>> genes = ProteinComplexAnalyzer.getUniprotGeneMapping()
							.mapUniprotACCToGeneByType(acc);
					if (!genes.isEmpty()) {

						Set<String> genesTMP = null;
						if (genes.containsKey(UniprotGeneMapping.GENE_NAME)) {
							genesTMP = genes.get(UniprotGeneMapping.GENE_NAME);
						} else {
							genesTMP = genes.values().iterator().next();
						}
						final StringBuilder sb = new StringBuilder();
						for (final String geneTMP : genesTMP) {
							if (!"".equals(sb.toString())) {
								sb.append(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
							}
							sb.append(geneTMP);
						}
						gene = sb.toString();
					}
					final ProteinComponent component = new ProteinComponent(acc, gene);
					complex.addComponent(component);
				}
			} else {
				// this is because it is a list of genes
				for (final String gene : components) {
					String acc = null;
					final Set<String> accs = ProteinComplexAnalyzer.getUniprotGeneMapping().mapGeneToUniprotACC(gene);
					if (!accs.isEmpty()) {
						final StringBuilder sb = new StringBuilder();
						for (final String acc2 : accs) {
							if (!"".equals(sb.toString())) {
								sb.append(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
							}
							sb.append(acc2);
						}
						acc = sb.toString();
					}
					final ProteinComponent component = new ProteinComponent(acc, gene);
					complex.addComponent(component);
				}
			}
		}
		// iterate over the complexes and reduce the accession to the Swissprot one if
		// possible
		final Set<String> toQuery = new THashSet<String>();
		for (final ProteinComplex proteinComplex : ret) {
			for (final ProteinComponent component : proteinComplex.getComponentList()) {
				if (component.getAcc().contains(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR)) {
					final String[] split = component.getAcc().split(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
					for (final String acc : split) {
						toQuery.add(acc);
					}
				}
			}
		}
		int numProteinsSimplified = 0;
		final Map<String, Entry> annotatedProteins = ProteinComplexAnalyzer.getUPLR().getAnnotatedProteins(null,
				toQuery);
		for (final ProteinComplex proteinComplex : ret) {
			for (final ProteinComponent component : proteinComplex.getComponentList()) {
				if (component.getAcc().contains(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR)) {
					final String[] split = component.getAcc().split(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
					final Set<String> accs = new THashSet<String>();
					for (final String acc : split) {
						accs.add(acc);
					}
					final Set<String> accsSimplified = getReviewed(accs, annotatedProteins);
					final StringBuilder accStringBuilder = new StringBuilder();
					for (final String acc : accsSimplified) {
						if (!"".equals(accStringBuilder.toString())) {
							accStringBuilder.append(ProteinComplexAnalyzer.AMBIGUOUS_SEPARATOR);
						}
						accStringBuilder.append(acc);
					}
					if (accsSimplified.size() < accs.size()) {
						numProteinsSimplified++;
					}
					component.setAcc(accStringBuilder.toString());
				}
			}
		}
		log.info(numProteinsSimplified
				+ " proteins in which we selected Swissprot ACC among all mapped ones from their gene name");

		showError(ret.size() + " protein complexes loaded", "Protein complexes loaded");
		return ret;
	}

	private List<ProteinComplex> loadProteinComplexesFromEPICFolder(File epicResultsFolder) throws IOException {
		final List<ProteinComplex> ret = EpicResultComparator.readPredictedProteinComplexes(epicResultsFolder);

		showError(ret.size() + " protein complexes loaded", "Protein complexes loaded");
		return ret;
	}

	private Set<String> getReviewed(Set<String> accs2, Map<String, Entry> annotatedProteins) {
		final Set<String> ret = new THashSet<String>();
		for (final String acc : accs2) {
			if (annotatedProteins.containsKey(acc)) {
				if (UniprotEntryUtil.isSwissProt(annotatedProteins.get(acc))) {
					ret.add(acc);
				}
			}
		}
		if (ret.isEmpty()) {
			return accs2;
		}
		return ret;
	}

	protected void addProfileData() {
		try {
			final JFileChooser fileChooser = new JFileChooser(currentDirectoryProfiles);
			fileChooser.showOpenDialog(this);
			final File selectedFile = fileChooser.getSelectedFile();
			currentDirectoryProfiles = selectedFile.getParentFile();
			experiment = loadProfile(selectedFile);
			showError("Loaded experiment with " + experiment.getFractions().size() + "  fractions",
					"Experiment loaded");
		} catch (final Exception e) {
			showError(e);
		}
	}

	private void showError(Exception e) {
		showError(e.getMessage(), "Error");

	}

	private void showError(String error, String title) {
		JOptionPane.showMessageDialog(this, error, title, JOptionPane.ERROR_MESSAGE);

	}

	private SeparationExperiment loadProfile(File selectedFile) throws IOException {
		final SeparationExperiment experiment = ProteinComplexAnalyzer
				.loadProjectSummaryFileNEW(FilenameUtils.getBaseName(selectedFile.getAbsolutePath()), selectedFile);
		return experiment;
	}

	public void showGraph(ProteinComplex selectedComplex, DataType dataType, boolean separated, Integer smoothWindow) {

		final String title = selectedComplex.toString(ProteinComplex.useGeneToPrint, "-");

		updateNumberOfProteinsOfSelectedComplex(selectedComplex.getComponents().size());

		graphContainerPanel.removeAll();
		int numRow = 0;
		if (!separated) {
			// spc
			final JFreeChart chart = getChart(selectedComplex.getComponentListName(ProteinComplex.useGeneToPrint),
					title.toString(), dataType, smoothWindow);
			final ChartPanel graphPanel = new ChartPanel(chart);
			graphPanel.setFillZoomRectangle(false);
			graphPanel.setBorder(new TitledBorder(null, "Chart", TitledBorder.LEADING, TitledBorder.TOP, null, null));
			final GridBagConstraints c = new GridBagConstraints();
			c.gridx = 0;
			c.gridy = numRow++;
			c.weightx = 1;
			c.weighty = 1;
			c.fill = GridBagConstraints.BOTH;
			graphContainerPanel.add(graphPanel, c);
		} else {
			for (final ProteinComponent component : selectedComplex.getComponentList()) {
				final List<String> accList = new ArrayList<String>();
				if (ProteinComplex.useGeneToPrint) {
					accList.add(component.getGene());
				} else {
					accList.add(component.getAcc());
				}
				final JFreeChart chart = getChart(accList, accList.get(0), dataType, smoothWindow);
				// only show title if the chart is just one
				chart.getTitle().setVisible(selectedComplex.getComponents().size() == 1);
				final ChartPanel graphPanel = new ChartPanel(chart);
				graphPanel.setFillZoomRectangle(false);
				graphPanel
						.setBorder(new TitledBorder(null, "Chart", TitledBorder.LEADING, TitledBorder.TOP, null, null));
				final GridBagConstraints c2 = new GridBagConstraints();
				c2.gridx = 0;
				c2.gridy = numRow++;
				c2.weightx = 1;
				c2.weighty = 1;
				c2.fill = GridBagConstraints.BOTH;
				graphContainerPanel.add(graphPanel, c2);
			}
		}
		setVisible(true);
	}

	private void updateNumberOfProteinsOfSelectedComplex(int size) {
		proteinsInSelectedClusterLabel.setText(size + " proteins");
	}

	/**
	 * 
	 * @param selectedComplex
	 * @param title
	 * @param smoothWindow
	 * @param showSPC         if true is SPC if false is NSAF
	 * @return
	 */
	public JFreeChart getChart(List<String> accs, String title, DataType dataType, Integer smoothWindow) {

		final NumberAxis numberaxis_x = new NumberAxis("Fraction");
		final String yLabel = dataType.name();
		final NumberAxis numberaxis_y = new NumberAxis(yLabel);

		final XYPlot xyplot = new XYPlot();
		xyplot.setDomainAxis(numberaxis_x);
		xyplot.setRangeAxis(numberaxis_y);
		// if (false) {
		// // this is to see vertical lines in each point
		// final XYBarRenderer xybarrenderer = new XYBarRenderer();
		// xybarrenderer.setShadowVisible(false);
		// xybarrenderer.setMargin(0.95);
		// xyplot.setRenderer(xybarrenderer);
		// } else {
		xyplot.setRenderer(new XYLineAndShapeRenderer());
		// }
		xyplot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
		int numDataSet = 1;
		for (final String acc : accs) {

			if (dataType == DataType.SPC) {
				// experimental data
				final IntervalXYDataset intervalxydataset = createDatasetFromSPCProfile(
						getSpcProfile(acc, smoothWindow), acc);
				xyplot.setDataset(numDataSet++, intervalxydataset);

				final NumberAxis axis = (NumberAxis) xyplot.getRangeAxis();
				// axis.setTickUnit(new NumberTickUnit(20));
				axis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

				// final XYBarRenderer xybarrenderer = new XYBarRenderer();
				// xybarrenderer.setShadowVisible(false);
				// xybarrenderer.setMargin(0.95);
				// xyplot.setRenderer(numDataSet - 1, xybarrenderer);
				xyplot.setRenderer(numDataSet - 1, new XYLineAndShapeRenderer());
			} else if (dataType == DataType.NSAF) {
				final IntervalXYDataset intervalxydataset = createDatasetFromNSAFProfile(
						getNSAFProfile(acc, smoothWindow), acc);
				xyplot.setDataset(numDataSet++, intervalxydataset);
				// final XYBarRenderer xybarrenderer = new XYBarRenderer();
				// xybarrenderer.setShadowVisible(false);
				// xybarrenderer.setMargin(0.95);
				// xyplot.setRenderer(numDataSet - 1, xybarrenderer);
				xyplot.setRenderer(numDataSet - 1, new XYLineAndShapeRenderer());
			}

		}

		final JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT, xyplot, true);
		final LegendTitle legend = chart.getLegend();
		legend.setPosition(RectangleEdge.RIGHT);

		return chart;
	}

	private List<double[]> getSpcProfile(String acc, Integer smoothWindow) {
		TDoubleArrayList spcs = new TDoubleArrayList();
		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Set<Protein> protein = fraction.getProteinByKey(acc);
			if (protein != null && protein.size() == 1) {
				spcs.add(protein.iterator().next().getSpc());
			} else {
				spcs.add(0);
			}
		}

		if (smoothWindow != null && smoothWindow > 0) {
			spcs = smooth(spcs, smoothWindow);
		}
		return PComplexUtil.transformExperimentalDataForModelling(spcs);
	}

	private TDoubleArrayList smooth(TDoubleArrayList spcs, Integer smoothWindow) {

		for (int index = 0; index < spcs.size(); index++) {
			final TDoubleArrayList toAverage = new TDoubleArrayList();
			for (int offset = -smoothWindow; offset <= smoothWindow; offset++) {
				if (index + offset >= 0 && index + offset < spcs.size()) {
					toAverage.add(spcs.get(index + offset));
				}
			}
			final double mean = Maths.mean(toAverage);
			spcs.set(index, mean);
		}
		return spcs;
	}

	private List<double[]> getNSAFProfile(String acc, Integer smoothWindow) {
		TDoubleArrayList nsafs = new TDoubleArrayList();
		for (final Fraction fraction : experiment.getSortedFractions()) {
			final Protein protein = fraction.getProteinByAcc(acc);
			if (protein != null) {
				nsafs.add(protein.getNSAF());
			} else {
				nsafs.add(0.0);
			}
		}
		if (smoothWindow != null && smoothWindow > 0) {
			nsafs = smooth(nsafs, smoothWindow);
		}
		return PComplexUtil.transformExperimentalDataForModelling(nsafs);
	}

	private IntervalXYDataset createDatasetFromSPCProfile(List<double[]> profile, String acc) {
		// log.info("Histogram values are:");
		// StringBuilder sb = new StringBuilder();
		// for (double[] values : experimentalHistogramData) {
		// sb.append("(" + values[0] + "," + values[1] + "), ");
		// }
		// // log.info(sb.toString());

		final XYSeries xyseries = new XYSeries(acc);
		for (int i = 0; i < profile.size(); i++) {
			final double x = profile.get(i)[0];
			final double y = profile.get(i)[1];
			xyseries.add(x, y);
		}
		return new XYSeriesCollection(xyseries);
	}

	private IntervalXYDataset createDatasetFromNSAFProfile(List<double[]> profile, String acc) {
		// log.info("Histogram values are:");
		// StringBuilder sb = new StringBuilder();
		// for (double[] values : experimentalHistogramData) {
		// sb.append("(" + values[0] + "," + values[1] + "), ");
		// }
		// // log.info(sb.toString());

		final XYSeries xyseries = new XYSeries(acc);
		for (int i = 0; i < profile.size(); i++) {
			final double x = profile.get(i)[0];
			final double y = profile.get(i)[1];
			xyseries.add(x, y);
		}
		return new XYSeriesCollection(xyseries);
	}

}
