package de.ovgu.mpa.validator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class FASTAWriter {

    private BufferedWriter bw;

    public FASTAWriter(String path) {
        try {
            bw = new BufferedWriter(new FileWriter(path));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void write(ArrayList<Protein> pL) throws IOException{
        for (Protein prot : pL) {
            bw.write(">" + prot.getDescription() + "\n");
            bw.write(prot.getSequence() + "\n");
        }
        bw.flush();
        bw.close();
    }
}