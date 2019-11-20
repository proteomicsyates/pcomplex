package edu.scripps.yates.pcomplex.gui;

import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;

import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.pcomplex.model.ProteinComponent;

public class ProteinComplexTable extends JTable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -2353166180758620517L;
	private static Logger log = Logger.getLogger("log4j.logger.org.proteored");
	private final static String[] headers = { "Acc", "Gene", "Name", "detected" };
	private static final int[] columnWidths = { 10, 10, 40, 4 };

	public ProteinComplexTable() {
		initTableColumns();
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

	private void addNewRow(List<String> lineStringList) {
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

	public void addNewComponent(ProteinComponent component, boolean detected) {
		final List<String> list = new ArrayList<String>();
		list.add(component.getAcc());
		list.add(component.getGene());
		list.add(component.getProteinName());
		list.add(String.valueOf(detected));
		addNewRow(list);
		SwingUtilities.invokeLater(new Runnable() {

			@Override
			public void run() {
				repaint();
			}
		});
	}

	public void setProteinComplex(ProteinComplex referenceComplex, ProteinComplex detectedComplex) {
		if (referenceComplex != null && detectedComplex != null) {
			clearData();
			for (final ProteinComponent component : referenceComplex.getComponentList()) {
				final boolean detected = detectedComplex.getComponents().containsKey(component.getKey());
				addNewComponent(component, detected);
			}
		}
	}

}
