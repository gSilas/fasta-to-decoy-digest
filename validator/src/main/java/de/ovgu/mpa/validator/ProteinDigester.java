package de.ovgu.mpa.validator;

import java.util.ArrayList;
import java.util.LinkedList;

public final class ProteinDigester {
    private static StringBuilder sb = new StringBuilder();

    public static LinkedList<String> digestProtein(String proteinSequence) {

		if(proteinSequence.contains("X") && ValidatorConfig.excludeX) {
			return new LinkedList<>();
		}

        LinkedList<String> peptides = new LinkedList<>();
		// simple digest, no missed cleavages, no length filtering
		ArrayList<String> peptidesMCZero = new ArrayList<String>(); 
		sb.setLength(0);
		for (char c : proteinSequence.toCharArray()) {
			sb.append(c);
			if (c == 'K' || c == 'R') {
				// change here for Peptide length
				peptidesMCZero.add(sb.toString());
				sb.setLength(0);
			}
		}
		// last peptide
		if (sb.length() != 0) {
			peptidesMCZero.add(sb.toString());	
		}
		// number of peptides minus 1 == max missed cleavages
		int missedCleavage =  peptidesMCZero.size() - 1;
		// consider missed cleavage (all of them)
		sb.setLength(0);
		for (int currentIndex = 0; currentIndex < peptidesMCZero.size(); currentIndex++) {
			for (int mc = 0; mc <= missedCleavage; mc++) {
				if (currentIndex + mc < peptidesMCZero.size()) {
					sb.append(peptidesMCZero.get(currentIndex + mc));
					if (sb.length() <  ValidatorConfig.MAXIMUM_PEP_LENGTH && sb.length() > ValidatorConfig.MINIMUM_PEP_LENGTH) {
						peptides.add(sb.toString());
					}
				}
			}
			sb.setLength(0);
        }
		return peptides;
    }
}