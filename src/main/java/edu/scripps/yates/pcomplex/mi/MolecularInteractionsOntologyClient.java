package edu.scripps.yates.pcomplex.mi;

import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.ols.web.service.client.OLSClient;
import uk.ac.ebi.pride.utilities.ols.web.service.config.OLSWsConfig;
import uk.ac.ebi.pride.utilities.ols.web.service.model.Identifier;
import uk.ac.ebi.pride.utilities.ols.web.service.model.Term;

public class MolecularInteractionsOntologyClient {
	private final static Logger log = Logger.getLogger(MolecularInteractionsOntologyClient.class);
	private static final OLSClient olsClient = new OLSClient(new OLSWsConfig());
	private static final String BIOCHEMICAL_OBO_ID = "MI:0401";

	public static boolean containsBiochemicalAsParent(String termID) {
		return containsAsParent(termID, BIOCHEMICAL_OBO_ID, true);
	}

	public static boolean containsAsParent(String termID, String parentID, boolean retryiffails) {
		try {
			final List<Term> parent = olsClient.getTermParents(new Identifier(termID, Identifier.IdentifierType.OBO),
					"mi", 100);
			if (parent.isEmpty()) {
				throw new IllegalArgumentException("term " + termID + " has no parents?? Maybe the id is wrong");
			}
			for (final Term term : parent) {
				if (term.getOboId().getIdentifier().equals(parentID)) {
					return true;
				}
			}

		} catch (final Exception e) {
			e.printStackTrace();
			log.error("Error getting term " + termID);
			if (retryiffails) {
				return containsAsParent(termID, parentID, false);
			}
		}
		return false;
	}

}
