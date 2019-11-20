package edu.scripps.yates.pcomplex.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.util.UniprotEntryEBIUtil;
import edu.scripps.yates.pcomplex.model.Protein;
import edu.scripps.yates.pcomplex.model.ProteinComplex;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

public class PComplexUtil {
	/**
	 * Gets a set of key strings that are the accessions and gene names of a
	 * collection of proteins
	 * 
	 * @param proteins
	 * @return
	 */
	public static Set<String> getKeysFromProteins(Collection<Protein> proteins) {
		final Set<String> keys = new THashSet<String>();
		for (final Protein protein : proteins) {
			keys.add(protein.getAcc());
			keys.add(protein.getGene());
			keys.addAll(protein.getOthers());
		}
		return keys;
	}

	/**
	 * Gets a set of key strings that are the keys of the proteins of the
	 * collection of {@link ProteinComplex}s
	 * 
	 * @param complexes
	 * @return
	 */
	public static Set<String> getKeysFromProteinComplexes(Collection<ProteinComplex> complexes) {
		final Set<String> keys = new THashSet<String>();
		for (final ProteinComplex complex : complexes) {
			keys.addAll(complex.getComponents().keySet());
		}
		return keys;
	}

	public static Map<String, Set<Protein>> getKeyMapFromProteins(Collection<Protein> proteins) {
		final Map<String, Set<Protein>> keys = new THashMap<String, Set<Protein>>();
		for (final Protein protein : proteins) {
			final Set<Protein> set = new THashSet<Protein>();
			set.add(protein);
			keys.put(protein.getAcc(), set);
			if (keys.containsKey(protein.getGene())) {
				keys.get(protein.getGene()).add(protein);
			} else {
				final Set<Protein> set2 = new THashSet<Protein>();
				set2.add(protein);
				keys.put(protein.getGene(), set2);
			}
			for (final String otherKey : protein.getOthers()) {
				if (keys.containsKey(otherKey)) {
					keys.get(otherKey).add(protein);
				} else {
					final Set<Protein> set2 = new THashSet<Protein>();
					set2.add(protein);
					keys.put(otherKey, set2);
				}
			}

		}
		return keys;
	}

	public static boolean hasTransmembraneRegion(String acc, UniprotProteinLocalRetriever uplr) {
		final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
		if (annotatedProtein.containsKey(acc)) {
			final Entry entry = annotatedProtein.get(acc);
			final List<edu.scripps.yates.utilities.annotations.uniprot.xml.FeatureType> features = UniprotEntryEBIUtil
					.getFeatures(entry, FeatureType.TRANSMEM);
			if (!features.isEmpty()) {
				return true;
			}
		}
		return false;
	}

	public static List<double[]> transformExperimentalDataForModelling(TDoubleList profile) {
		final List<double[]> ret = new ArrayList<double[]>();
		int x = 1;
		for (final double y : profile.toArray()) {
			final double[] point = new double[2];
			point[0] = x;
			point[1] = y;
			ret.add(point);
			x++;
		}
		return ret;
	}

	public static List<int[]> transformExperimentalDataForModelling(TIntList profile) {
		final List<int[]> ret = new ArrayList<int[]>();
		int x = 1;
		for (final int y : profile.toArray()) {
			final int[] point = new int[2];
			point[0] = x;
			point[1] = y;
			ret.add(point);
			x++;
		}
		return ret;
	}

}
