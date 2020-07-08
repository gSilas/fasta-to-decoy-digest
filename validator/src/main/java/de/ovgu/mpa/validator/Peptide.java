package de.ovgu.mpa.validator;

import java.util.Comparator;

public class Peptide {
	public Peptide(String s, double m) {
		this.sequence = s;
		this.mass = m;
	}
	public String sequence;
	public double mass;
	
	public static Comparator<Peptide> getComparator() {
		return (new Comparator<Peptide>() {
		    @Override
		    public int compare(Peptide c1, Peptide c2) {
		        return Double.compare(c1.mass, c2.mass);
		    }
		});
	}
	
}