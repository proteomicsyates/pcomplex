package edu.scripps.yates.pcomplex.gprofiler;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.log4j.Logger;

public class GProfiler {
	private final static Logger log = Logger.getLogger(GProfiler.class);

	public static GProfilerResult readGProfilerResults(File gProfilerResult) throws IOException {
		final GProfilerResult ret = new GProfilerResult();
		final List<String> lines = Files.readAllLines(gProfilerResult.toPath());

		for (final String line : lines) {

			if (line.startsWith("\"GO")) {
				final GOTerm go = new GOTerm(line.split("\",\""));
				ret.addGOTerm(go);
			}
		}
		log.info(ret.getNumGOTerms() + " GOTerms read from file " + gProfilerResult.getAbsolutePath());
		return ret;
	}
}
