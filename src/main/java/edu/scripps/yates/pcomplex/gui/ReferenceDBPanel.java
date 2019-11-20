package edu.scripps.yates.pcomplex.gui;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.border.TitledBorder;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;

public class ReferenceDBPanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = -6061913891748092349L;
	private static final double MIN_OVERLAP = 0.3;

	private final List<ProteinComplexDB> dbs = new ArrayList<ProteinComplexDB>();

	private final ProteinComplexOverlapTable overlapTable;
	private final ProteinComplexTable complexesTable;
	private final JPanel headerPanel;
	private ProteinComplex selectedDetectedComplex;

	public ReferenceDBPanel() {
		setLayout(new BorderLayout(0, 0));
		// setBorder(new TitledBorder(null, db.getName(), TitledBorder.LEADING,
		// TitledBorder.TOP, null, null));
		headerPanel = new JPanel();
		add(headerPanel, BorderLayout.NORTH);
		headerPanel.setLayout(new GridBagLayout());

		final JPanel panel_3 = new JPanel();
		add(panel_3);
		panel_3.setLayout(new GridLayout(0, 1, 0, 0));

		final JPanel panel_1 = new JPanel();
		panel_3.add(panel_1);
		panel_1.setBorder(
				new TitledBorder(null, "Overlapped complexes", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_1.setLayout(new BorderLayout(0, 0));

		complexesTable = new ProteinComplexTable();
		overlapTable = new ProteinComplexOverlapTable(this);
		overlapTable.setAutoResizeMode(JTable.AUTO_RESIZE_NEXT_COLUMN);
		final JScrollPane scroll = new JScrollPane(overlapTable);
		panel_1.add(scroll);

		final JPanel panel_2 = new JPanel();
		panel_3.add(panel_2);
		panel_2.setBorder(new TitledBorder(null, "Selected reference complex", TitledBorder.LEADING, TitledBorder.TOP,
				null, null));
		panel_2.setLayout(new BorderLayout(0, 0));

		final JScrollPane scrollPane = new JScrollPane(complexesTable);
		panel_2.add(scrollPane, BorderLayout.CENTER);

	}

	public void addDB(ProteinComplexDB db) {
		dbs.add(db);
		final GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = dbs.size() - 1;
		final JLabel lblName = new JLabel("Name:");
		lblName.setHorizontalAlignment(SwingConstants.RIGHT);
		headerPanel.add(lblName, c);

		c.gridx = 1;
		final JLabel dbNameLabel = new JLabel(db.getName());
		headerPanel.add(dbNameLabel, c);

		c.gridx = 2;
		final JLabel lblNewLabel = new JLabel("# complexes:");
		lblNewLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		headerPanel.add(lblNewLabel, c);

		c.gridx = 3;
		final JLabel numComplexesLabel = new JLabel(String.valueOf(db.getProteinComplexes().size()));
		headerPanel.add(numComplexesLabel, c);

	}

	public void loadOverlapsWithComplex(ProteinComplex selectedDetectedComplex) {
		overlapTable.clearData();
		this.selectedDetectedComplex = selectedDetectedComplex;
		for (final ProteinComplexDB db : dbs) {
			for (final ProteinComplex referenceComplex : db.getProteinComplexes()) {
				final double overlap = ClusterEvaluation.getOverlap(referenceComplex, selectedDetectedComplex);
				if (overlap > MIN_OVERLAP) {
					overlapTable.addNewComplex(db, referenceComplex, selectedDetectedComplex, overlap);
				}
			}
		}
	}

	public ProteinComplexTable getcomplexesTable() {
		return complexesTable;
	}

	public ProteinComplex getSelectedDetectedComplex() {
		return selectedDetectedComplex;

	}

	public ProteinComplex getComplexFromDB(String databaseName, String complexID) {
		for (final ProteinComplexDB db : dbs) {
			if (db.getName().equals(databaseName)) {
				return db.getProteinComplexesByID(complexID);
			}
		}
		return null;
	}
}
