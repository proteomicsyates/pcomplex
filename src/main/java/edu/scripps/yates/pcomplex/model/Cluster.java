package edu.scripps.yates.pcomplex.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import gnu.trove.set.hash.THashSet;
import uk.ac.rhul.cs.cl1.ValuedNodeSet;

/**
 * Set of {@link Protein} or {@link ProteinComponent} that are clustered
 * together because they interact
 * 
 * @author salvador
 *
 */
public class Cluster {

	private final ValuedNodeSet nodeSet;
	private ArrayList<ProteinComponent> proteinComponents;

	public Cluster(ValuedNodeSet nodeSet) {
		this.nodeSet = nodeSet;
	}

	public String getProteinNames(String separator) {
		final Set<String> set = new THashSet<String>();
		final StringBuilder sb = new StringBuilder();
		for (final String member : nodeSet.getMemberNames()) {
			if (set.contains(member)) {
				continue;
			}
			set.add(member);
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(member);
		}
		return sb.toString();
	}

	public double getSignificance() {
		return nodeSet.getSignificance();
	}

	public double getDensity() {
		return nodeSet.getDensity();
	}

	public List<ProteinComponent> getProteinComponents() throws IOException {
		if (proteinComponents == null) {
			proteinComponents = new ArrayList<ProteinComponent>();
			final Set<String> set = new THashSet<String>();
			for (final String member : nodeSet.getMemberNames()) {
				if (set.contains(member)) {
					continue;
				}
				proteinComponents.add(new ProteinComponent(member, null));
			}
		}
		return proteinComponents;
	}
}
