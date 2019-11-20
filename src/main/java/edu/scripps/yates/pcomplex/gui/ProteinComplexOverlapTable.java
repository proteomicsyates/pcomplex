package edu.scripps.yates.pcomplex.gui;

import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.db.ProteinComplexDB;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.util.ClusterEvaluation;

public class ProteinComplexOverlapTable extends JTable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -2353166180758620517L;
	private static Logger log = Logger.getLogger("log4j.logger.org.proteored");
	private final static String DB = "DB";
	private final static String ID = "ID";
	private static final String[] headers = { DB, ID, "Name", "# components", "overlap score", "# common components" };
	private static final int[] columnWidths = { 20, 10, 40, 10, 10, 10 };

	public ProteinComplexOverlapTable(ReferenceDBPanel referencePanel) {
		setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		setRowSelectionAllowed(true);
		setColumnSelectionAllowed(false);
		initTableColumns();
		final ListSelectionModel cellSelectionModel = getSelectionModel();
		cellSelectionModel.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

		cellSelectionModel.addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e) {
				final String selectedData = null;

				final int[] selectedRow = getSelectedRows();
				final int[] selectedColumns = getSelectedColumns();
				referencePanel.getcomplexesTable().clearData();
				for (int i = 0; i < selectedRow.length; i++) {
					log.info("Selected row: " + selectedRow[i]);
					final String complexID = (String) getValueAt(selectedRow[i], Arrays.asList(headers).indexOf(ID));
					final String databaseName = (String) getValueAt(selectedRow[i], Arrays.asList(headers).indexOf(DB));
					final ProteinComplex referenceComplex = referencePanel.getComplexFromDB(databaseName, complexID);
					referencePanel.getcomplexesTable().setProteinComplex(referenceComplex,
							referencePanel.getSelectedDetectedComplex());

				}

			}

		});
	}

	@Override
	protected JTableHeader createDefaultTableHeader() {
		return new JTableHeader(columnModel) {
			@Override
			public String getToolTipText(MouseEvent e) {
				final java.awt.Point p = e.getPoint();
				final int index = columnModel.getColumnIndexAtX(p.x);
				// int realIndex =
				// columnModel.getColumn(index).getModelIndex();
				final String columnName = (String) columnModel.getColumn(index).getHeaderValue();
				final String tip = getToolTipTextForColumn(columnName);
				// log.info("Tip = " + tip);
				if (tip != null)
					return tip;
				else
					return super.getToolTipText(e);
			}
		};
	}

	private String getToolTipTextForColumn(String columnName) {
		return columnName;
	}

	public void clearData() {

		final TableModel model = getModel();
		if (model instanceof DefaultTableModel) {
			((DefaultTableModel) model).setRowCount(0);
			// ((DefaultTableModel) model).setColumnCount(0);
		}

	}

	private void addColumnsInTable(List<String> columnsStringList) {
		final DefaultTableModel defaultModel = getTableModel();
		log.info("Adding colums " + columnsStringList.size() + " columns");
		if (columnsStringList != null) {

			for (final String columnName : columnsStringList) {
				defaultModel.addColumn(columnName);
			}
			log.info("Added " + getColumnCount() + " colums");
			for (int i = 0; i < getColumnCount(); i++) {
				final TableColumn column = getColumnModel().getColumn(i);

				column.setPreferredWidth(columnWidths[i]);

				column.setResizable(true);
			}
		}
	}

	private DefaultTableModel getTableModel() {

		TableModel model = getModel();
		if (model == null) {
			model = new DefaultTableModel();
		}
		final DefaultTableModel defaultModel = (DefaultTableModel) model;
		return defaultModel;
	}

	public void addNewRow(List<String> lineStringList) {
		final DefaultTableModel defaultModel = getTableModel();
		try {
			defaultModel.addRow(lineStringList.toArray());
		} catch (final IndexOutOfBoundsException e) {
			e.printStackTrace();
			log.error(lineStringList);
		}
		// this.table.repaint();

	}

	private void initTableColumns() {
		addColumnsInTable(Arrays.asList(headers));
	}

	public void addNewComplex(ProteinComplexDB db, ProteinComplex referenceComplex, ProteinComplex discoveredComplex,
			double overlap) {
		log.info("Adding new row to table with complex: " + referenceComplex.getComponentsKey());
		final List<String> list = new ArrayList<String>();
		list.add(db.getName());
		list.add(referenceComplex.getId());
		list.add(referenceComplex.getName());
		list.add(String.valueOf(referenceComplex.getComponentList().size()));
		list.add(String.valueOf(overlap));
		list.add(String.valueOf(ClusterEvaluation.getIntersection(referenceComplex, discoveredComplex).size()));
		addNewRow(list);
		SwingUtilities.invokeLater(new Runnable() {

			@Override
			public void run() {
				repaint();
			}
		});
	}
}
