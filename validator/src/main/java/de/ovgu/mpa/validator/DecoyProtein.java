package de.ovgu.mpa.validator;

import java.util.LinkedList;

public class DecoyProtein extends Protein {
	
	public DecoyProtein(FastaProtein prot) {
        this.peptides = new LinkedList<>();
        StringBuilder sb = new StringBuilder(prot.getSequence());
        sb.reverse();
        this.description = "DECOY_" + prot.getDescription();
        this.sequence = sb.toString();
        this.peptides = ProteinDigester.digestProtein(this.sequence);
	}

}