package de.ovgu.mpa.validator.fasta;

import java.util.LinkedList;

import de.ovgu.mpa.validator.Protein;
import de.ovgu.mpa.validator.ProteinDigester;

public class FastaProtein extends Protein {
	
	public FastaProtein() {
		this.peptides = new LinkedList<>();
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
		this.peptides = ProteinDigester.digestProtein(this.sequence);
	}

}
