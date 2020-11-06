package de.ovgu.mpa.validator.mascot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import de.ovgu.mpa.validator.AminoAcid;
import de.ovgu.mpa.validator.Constants;
import de.ovgu.mpa.validator.ValidatorConfig;

public class MascotCSVReader {

    BufferedReader csvReader;

    public MascotCSVReader(File csv) {
        try {
            this.csvReader = new BufferedReader(new FileReader(csv));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    public long getPeptides(String nonRedundantPepDB) throws IOException {

        BufferedWriter bwNR = new BufferedWriter(new FileWriter(new File(nonRedundantPepDB)));

        long resultCount = 0;

        System.out.println("Processing mascot csv!");

        String line = this.csvReader.readLine();

        System.out.println("Searching queries");

        while (!line.contains("\"query_number\",\"moverz\",\"charge\",\"intensity\",\"TotalIonsIntensity\"")) {
            line = this.csvReader.readLine();
        }

        System.out.println("Processing queries");

        line = this.csvReader.readLine();

        boolean textString = false;
        String currentString = "";
        int j = 0;
        String peptideIons1 = null;
        String peptideScore = null;
        String peptideSequence = null;

        double peptideMass = 0.0;
        double masswater = Constants.MASS_WATER;

        for (int i = 0; i < line.length(); i++) {
            if (line.charAt(i) == ',') {
                if (!textString) {
                    if (j == 6){
                        peptideIons1 = currentString;

                    } else if (j == 13){
                        peptideScore = currentString;

                    } else if (j == 15) {
                        peptideSequence = currentString;

                    }
                    currentString = "";
                    j++;
                }
                else if (textString) {
                    currentString += line.charAt(i);
                }
            }
            else if (line.charAt(i) == '\"') {
                textString = !textString;
            } else {
                currentString += line.charAt(i);
            }
        }

        if (peptideIons1 != null && peptideScore != null && peptideSequence != null) {
            if (!peptideIons1.isEmpty() && !peptideScore.isEmpty() && !peptideSequence.isEmpty()) {
                if (peptideSequence.length() > ValidatorConfig.MINIMUM_PEP_LENGTH && peptideSequence.length() < ValidatorConfig.MAXIMUM_PEP_LENGTH) {
                    resultCount++;
                    peptideMass = 0.0;
                    for (char c : peptideSequence.toCharArray()) {
                        AminoAcid aa = AminoAcid.fromOneLetter(c);
                        peptideMass += aa.getMass();
                    }
                    peptideMass -= masswater * (peptideSequence.length() - 1);
    
                    bwNR.write(peptideSequence + ";" + Double.toString(peptideMass) + ";;" + peptideScore + ";" + peptideIons1 + "\n");    
                }
            }
        }

        peptideIons1 = null;
        peptideScore = null;
        peptideSequence = null;

        while (true) {
            line = this.csvReader.readLine();

            if (line == null ) {
                break;
            }

            textString = false;
            currentString = "";
            j = 0;
            for (int i = 0; i < line.length(); i++) {
                if (line.charAt(i) == ',') {
                    if (!textString) {
                        if (j == 6){
                            peptideIons1 = currentString;

                        } else if (j == 13){
                            peptideScore = currentString;

                        } else if (j == 15) {
                            peptideSequence = currentString;

                        }
                        currentString = "";
                        j++;
                    }
                    else if (textString) {
                        currentString += line.charAt(i);
                    }
                }
                else if (line.charAt(i) == '\"') {
                    textString = !textString;
                } else {
                    currentString += line.charAt(i);
                }
            }

            if (peptideIons1 != null && peptideScore != null && peptideSequence != null) {
                if (!peptideIons1.isEmpty() && !peptideScore.isEmpty() && !peptideSequence.isEmpty()) {
                    if (peptideSequence.length() > ValidatorConfig.MINIMUM_PEP_LENGTH && peptideSequence.length() < ValidatorConfig.MAXIMUM_PEP_LENGTH) {
                        resultCount++;
                        peptideMass = 0.0;
                        for (char c : peptideSequence.toCharArray()) {
                            AminoAcid aa = AminoAcid.fromOneLetter(c);
                            peptideMass += aa.getMass();
                        }
                        peptideMass -= masswater * (peptideSequence.length() - 1);
        
                        bwNR.write(peptideSequence + ";" + Double.toString(peptideMass) + ";;" + peptideScore + ";" + peptideIons1 + "\n");    
                    }
                }
            }

            peptideIons1 = null;
            peptideScore = null;
            peptideSequence = null;
        }

        bwNR.flush();
        bwNR.close();
        
        return resultCount;
    }
    
}
