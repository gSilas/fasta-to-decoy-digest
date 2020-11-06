package de.ovgu.mpa.validator;

import java.util.LinkedList;

import de.ovgu.mpa.validator.fasta.FastaProtein;

public class DecoyProtein extends Protein {
	
	public DecoyProtein(FastaProtein fastaProtein) {
        this.peptides = new LinkedList<>();
        StringBuilder sb = new StringBuilder(fastaProtein.getSequence());
        sb.reverse();
        this.description = "DECOY_" + fastaProtein.getDescription();
        this.sequence = sb.toString();
        this.peptides = ProteinDigester.digestProtein(this.sequence);
	}

}