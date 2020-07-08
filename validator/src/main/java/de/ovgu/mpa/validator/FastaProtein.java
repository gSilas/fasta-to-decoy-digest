package de.ovgu.mpa.validator;

import java.util.LinkedList;

public class FastaProtein extends Protein{
	
	public FastaProtein() {
		this.peptides = new LinkedList<>();
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
		this.peptides = ProteinDigester.digestProtein(this.sequence);
	}

}
