package de.ovgu.mpa.validator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

public class DecoyGenerator {

    public static void generateReverseDecoyDatabase(Path fastaPath, Path decoyPath) {
        try {

            FASTAFileReader fr = new FASTAFileReader(fastaPath.toFile());
            fr.open();
            BufferedWriter bw = new BufferedWriter(new FileWriter(decoyPath.toFile()));
            while (fr.hasNext()) {
                FastaProtein fastaProtein = fr.next();
                DecoyProtein prot = new DecoyProtein(fastaProtein);
                bw.write(">" + prot.getDescription() + "\n");
                bw.write(prot.getSequence() + "\n");
            }

            fr.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void generateRandomDecoyDatabase(Path fastaPath, Path decoyPath) {
        // Requires aminoacid distribution
    }
}