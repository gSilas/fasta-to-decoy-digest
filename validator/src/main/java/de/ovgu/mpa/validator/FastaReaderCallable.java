package de.ovgu.mpa.validator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.Callable;

public class FastaReaderCallable implements Callable<Long> {

    String batchFiles;
    File fasta;
    long threadNumber;
    long startProtein;
    long endProtein;

    public FastaReaderCallable(String batchFiles, File fasta, long threadNumber, long startProtein, long endProtein) {
        this.batchFiles = batchFiles;
        this.fasta = fasta;
        this.threadNumber = threadNumber;

        this.startProtein = startProtein;
        this.endProtein = endProtein;
    }

    @Override
    public Long call() throws Exception {
        double masswater = Constants.MASS_WATER;
        
        long proteinCount = 0;
		double peptideMass;

		ArrayList<Peptide> databasePeptideList = new ArrayList<Peptide>();
		FASTAFileReader fr = new FASTAFileReader(fasta);
		fr.open();
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(batchFiles + Long.toString(threadNumber) + ".pep")));

        Protein fastaProtein = null;
        
        while (fr.hasNext()) {
			fastaProtein = fr.next();
            proteinCount++;
            if (proteinCount >= startProtein) break;
        }

        while ( proteinCount < endProtein && fr.hasNext()) {
            LinkedList<String> peptideList = fastaProtein.getPeptides();
            for (String peptideString : peptideList) {
                peptideMass = 0.0;
                for (char c : peptideString.toCharArray()) {
                    AminoAcid aa = AminoAcid.fromOneLetter(c);
                    peptideMass += aa.getMass();
                }
                
                peptideMass -= masswater * (peptideString.length() - 1);
                databasePeptideList.add(new Peptide(peptideString, peptideMass));
            }

            fastaProtein = fr.next();
            proteinCount++;
        }
        fr.close();
        
        // this step is crucial, it sorts the peptides by their mass enabling easy
        // merging later
        databasePeptideList.sort(Peptide.getComparator());
        for (Peptide peptide : databasePeptideList) {
            bw.write(peptide.sequence + ";" + Double.toString(peptide.mass) + "\n");
        }
        databasePeptideList.clear();
        bw.flush();
        bw.close();

		return proteinCount;
    }
    
}
