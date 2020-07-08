package de.ovgu.mpa.validator;

import java.util.LinkedList;

public class DecoyProtein extends Protein {
	
	public DecoyProtein(FastaProtein fp) {
        this.peptides = new LinkedList<>();
        StringBuilder sb = new StringBuilder(fp.getSequence());
        sb.reverse();
        this.description = "DECOY " + fp.getDescription();
        this.sequence = sb.toString();
        this.peptides = ProteinDigester.digestProtein(this.sequence);
	}

}