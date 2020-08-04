package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class PeptideReader {

    public PeptideReader() {

    }

    public ArrayList<Peptide> readPeptides(String pathToPeptides) throws IOException {
        BufferedReader brR = new BufferedReader(new FileReader(new File(pathToPeptides)));  
        ArrayList<Peptide> peptides = new ArrayList<Peptide>();
        String line = brR.readLine();
		while (line != null) {
			String[] split = line.split(";");
            String sequence = split[0];
            Double mass = Double.valueOf(split[1]);
            peptides.add(new Peptide(sequence, mass));
			line = brR.readLine();
		}
        brR.close();
        return peptides;
    }
}