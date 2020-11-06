package de.ovgu.mpa.validator.fasta;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import de.ovgu.mpa.validator.Protein;

public class FASTAWriter {

    private BufferedWriter bw;

    public FASTAWriter(String path) {
        try {
            bw = new BufferedWriter(new FileWriter(path));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void write(ArrayList<? extends Protein> pL) throws IOException{
        for (Protein prot : pL) {
            bw.write(">" + prot.getDescription() + "\n");
            bw.write(prot.getSequence() + "\n");
        }
        bw.flush();
        bw.close();
    }
}