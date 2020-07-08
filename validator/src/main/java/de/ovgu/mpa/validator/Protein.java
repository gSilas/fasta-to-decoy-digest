package de.ovgu.mpa.validator;

import java.util.LinkedList;

abstract public class Protein {

	protected String description;
    protected String sequence;

    protected LinkedList<String> peptides;

	public void setSequence(String sequence) {
		this.sequence = sequence;
    }
	
	public LinkedList<String> getPeptides() {
		return this.peptides;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}
	
	public String getSequence() {
		return sequence;
	}
}