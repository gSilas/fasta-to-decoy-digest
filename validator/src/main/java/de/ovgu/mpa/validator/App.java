package de.ovgu.mpa.validator;

import java.io.*;
import java.util.ArrayList;

public class App {
    public static void main(String[] args) {
        File file = new File(
                "/Users/dan/Downloads/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.fasta");
        /*
         * try { PeptideGenerator proteinDigester = new PeptideGenerator();
         * proteinDigester.createFiles("/Users/dan/results/", "/Users/dan/batches/",
         * file); } catch (IOException e) { e.printStackTrace(); }
         */

        FASTAFileReader fReader;
        try {
            fReader = new FASTAFileReader(file);
            fReader.open();
            ArrayList<Protein> pList = new ArrayList<Protein>();
            while(fReader.hasNext()) {
                FastaProtein prot = fReader.next();
                pList.add(new DecoyProtein(prot));
            }
            FASTAWriter fWriter = new FASTAWriter("/Users/dan/meme.fasta");
            fWriter.write(pList);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
