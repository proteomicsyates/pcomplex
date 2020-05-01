package edu.scripps.yates.pcomplex.gprofiler;

public enum GO_TYPE {
	CC, BP, MF;

	public static GO_TYPE parseString(String goTypeString) {
		if ("GO:MF".equals(goTypeString)) {
			return MF;
		} else if ("GO:BP".equals(goTypeString)) {
			return BP;
		} else if ("GO:CC".equals(goTypeString)) {
			return CC;
		}
		throw new IllegalArgumentException(goTypeString + " not recognized as a valid GO_Type");
	}
}
