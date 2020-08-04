package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

public class Statistics {

    public Statistics() {
    }

    public void comparePeptides(final String decoyPeptides, final String targetPeptides, String folder) throws IOException {
        final BufferedReader brD = new BufferedReader(new FileReader(new File(decoyPeptides)));
        final BufferedReader brT = new BufferedReader(new FileReader(new File(targetPeptides)));
        final double tolerance = 0.1;
        final double greatMatchTolerance = 0.5;
        int matchCounter = 0;
        int greatMatchCounter = 0;

        // <length, target proteins>
        int[] lengthTargetArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];

        // <length, decoy proteins>
        int[] lengthDecoyArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];

        // <length, bins>
        int[][] lengthBinMatrix = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH][10];

        // <length, great matches>
        // HashMap<Integer, Integer> lengthMap = new HashMap<>();
        int[] lengthMatchArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];

        for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++){
            lengthDecoyArray[i] = 0;
            lengthTargetArray[i] = 0;
            lengthMatchArray[i] = 0;
            for (int j = 0; j < 10; j++){
                lengthBinMatrix[i][j] = 0;
            }
        }

        String targetLine = brT.readLine();
        String[] targetSplit = targetLine.split(";");
        String targetSequence = targetSplit[0];
        Double targetMass = Double.valueOf(targetSplit[1]);

        lengthTargetArray[targetSequence.length()]++;

        String decoyLine = brD.readLine();
        String[] decoySplit = decoyLine.split(";");
        String decoySequence = decoySplit[0];
        Double decoyMass = Double.valueOf(decoySplit[1]);

        lengthDecoyArray[decoySequence.length()]++;

        final List<String[]> scoringDataLines = new ArrayList<>();
        scoringDataLines.add(new String[] { "Target", "Decoy", "TargetMass", "DecoyMass", "Cosine Similarity", "Euclidian Distance"});

        // <FragmentMatches, Occurences>
        HashMap<Integer, Integer> fragmentMatchMap = new HashMap<>();

        int[] cosineSimilarityBins = new int[10];
        for (int i = 0; i < cosineSimilarityBins.length; i++){
            cosineSimilarityBins[i] = 0;
        }

        boolean decoyChanged = false;
        boolean targetChanged = false;
        
        while (targetLine != null && decoyLine != null) {

            targetSplit = targetLine.split(";");
            targetSequence = targetSplit[0];
            targetMass = Double.valueOf(targetSplit[1]);

            if (targetChanged) lengthTargetArray[targetSequence.length()]++; targetChanged = false;

            decoySplit = decoyLine.split(";");
            decoySequence = decoySplit[0];
            decoyMass = Double.valueOf(decoySplit[1]);

            if (decoyChanged) lengthDecoyArray[targetSequence.length()]++; decoyChanged = false;
            
            if (targetMass <= decoyMass + tolerance && targetMass >= decoyMass - tolerance) {
                // MATCH
                HashSet<FragmentationCalculator.IonSeries> iontypes = new HashSet<FragmentationCalculator.IonSeries>();
                iontypes.add(FragmentationCalculator.IonSeries.B_ION_1);
                iontypes.add(FragmentationCalculator.IonSeries.Y_ION_1);
                iontypes.add(FragmentationCalculator.IonSeries.B_ION_2);
                iontypes.add(FragmentationCalculator.IonSeries.Y_ION_2);
                double[] targetIons = FragmentationCalculator.getFragmentIons(iontypes, targetSequence);
                double[] decoyIons = FragmentationCalculator.getFragmentIons(iontypes, decoySequence);
                
                ArrayList<Double> targetSpectrum = new ArrayList<Double>();
                ArrayList<Double> decoySpectrum = new ArrayList<Double>();

                int indexDecoy = 0;
                int indexTarget = 0;
                int fragmentMatch = 0;

                while(indexDecoy < decoyIons.length && indexTarget < targetIons.length) {
                    double targetIon = targetIons[indexTarget];
                    double decoyIon = decoyIons[indexDecoy];
                    
                    if (targetIon <= decoyIon + tolerance && targetIon >= decoyIon - tolerance){
                        targetSpectrum.add(1.);
                        decoySpectrum.add(1.);
                        fragmentMatch++;
                        indexTarget++;
                        indexDecoy++;
                    } else if (targetIon > decoyIon + tolerance) {
                        targetSpectrum.add(0.);
                        decoySpectrum.add(1.);
                        indexDecoy++;
                    } else if (targetIon < decoyIon - tolerance) {
                        targetSpectrum.add(1.);
                        decoySpectrum.add(0.);
                        indexTarget++;
                    }
                }

                while(indexDecoy < decoyIons.length){
                    targetSpectrum.add(0.);
                    decoySpectrum.add(1.);
                    indexDecoy++;
                }

                while(indexTarget < targetIons.length){
                    targetSpectrum.add(1.);
                    decoySpectrum.add(0.);
                    indexTarget++;
                }

                assert(decoySpectrum.size() == targetSpectrum.size());

                double cosineSimilarity = this.cosineSimilarity(this.toPrimitive(targetSpectrum), this.toPrimitive(decoySpectrum));
                
                if (cosineSimilarity > greatMatchTolerance) {
                    lengthMatchArray[targetSequence.length()]++;
                    greatMatchCounter++;
                }

                if (fragmentMatchMap.get(fragmentMatch) != null) {
                    fragmentMatchMap.put(fragmentMatch, fragmentMatchMap.get(fragmentMatch) + 1);
                } else {
                    fragmentMatchMap.put(fragmentMatch, 1);
                }

                int index = (int)(Math.round(cosineSimilarity * 10) % 10);
                if (Math.round(cosineSimilarity * 10) > 9.5) {
                    cosineSimilarityBins[9]++; 
                    lengthBinMatrix[targetSequence.length()][9]++;
                } else {
                    cosineSimilarityBins[index]++;
                    lengthBinMatrix[targetSequence.length()][index]++;
                }


                double euclidianDistance = this.ndistance(this.toPrimitive(targetSpectrum), this.toPrimitive(decoySpectrum));

                scoringDataLines.add(new String[] {targetSequence, decoySequence, Double.toString(targetMass), Double.toString(decoyMass), Double.toString(cosineSimilarity), Double.toString(euclidianDistance)});
                
                decoyLine = brD.readLine();
                decoyChanged = true;

                targetLine = brT.readLine();
                targetChanged = true;

                matchCounter++;

            }
            else if (targetMass > decoyMass + tolerance){
                decoyLine = brD.readLine();
                decoyChanged = true;

            }
            else if (targetMass < decoyMass - tolerance) {
                targetLine = brT.readLine();
                targetChanged = true;
            }
        }

        final List<String[]> lengthMatchLines = new ArrayList<>();
        final List<String[]> lengthCollectionLines = new ArrayList<>();

        lengthMatchLines.add(new String[] { "length", "matches" });
        lengthCollectionLines.add(new String[] { "length", "target", "peptides", "matches", "bins" });


        for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++){
            lengthMatchLines.add(new String[] {
                Integer.toString(i), 
                Integer.toString(lengthMatchArray[i])
            });

            String bString = "";
            for (int j = 0 ; j < 10; j++){
                bString += "," + lengthBinMatrix[i][j];
            }
            lengthCollectionLines.add(new String[] {
                Integer.toString(i), 
                Integer.toString(lengthTargetArray[i]), 
                Integer.toString(lengthDecoyArray[i]), 
                Integer.toString(lengthMatchArray[i]), 
                bString
            });
        }
        final List<String[]> fragmentMatchLines = new ArrayList<>();
        fragmentMatchLines.add(new String[] { "matches", "relative occurences" });

        for (final Map.Entry<Integer, Integer> entry : fragmentMatchMap.entrySet()) {
            fragmentMatchLines.add(new String[] { entry.getKey().toString(), Double.toString(100*((double) entry.getValue() / (double) matchCounter))});
        }

        final List<String[]> cosineSimilarityBinsLines = new ArrayList<>();
        cosineSimilarityBinsLines.add(new String[] { "last digit after comma", "occurences" });
        for (int i = 0; i < cosineSimilarityBins.length; i++) {
            cosineSimilarityBinsLines.add(new String[] {Integer.toString(i), Integer.toString(cosineSimilarityBins[i])});
        }

        final CSVWriter writer = new CSVWriter();
        try {
            writer.createCSV(scoringDataLines, folder + "/" + "scoring.csv");
            writer.createCSV(fragmentMatchLines, folder + "/" + "fragmentMatch.csv");
            writer.createCSV(lengthMatchLines, folder + "/" + "lenghtMatch.csv");
            writer.createCSV(cosineSimilarityBinsLines, folder + "/" + "cosineSimilarityBins.csv");
            writer.createCSV(lengthCollectionLines, folder + "/" + "lengthCollection.csv");
        } catch (final FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        System.out.println("Recognized " + matchCounter + " matches!");
        System.out.println("Recognized " + greatMatchCounter + " great matches!");

        brT.close();
        brD.close();
    }

    public double cosineSimilarity(double[] vectorA, double[] vectorB) {
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            dotProduct += vectorA[i] * vectorB[i];
            normA += Math.pow(vectorA[i], 2);
            normB += Math.pow(vectorB[i], 2);
        }   
        return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
    }

    public double ndistance(double[] a, double[] b) {
        double total = 0, diff;
        for (int i = 0; i < a.length; i++) {
            diff = b[i] - a[i];
            total += diff * diff;
        }
        return Math.sqrt(total);
    }

    public double[] padWithZero (double[] vector, int length) {
        int oldIndex = vector.length;
        vector = Arrays.copyOf(vector, length);
        Arrays.fill(vector, oldIndex, vector.length - 1, 0.0);
        return vector;
    }

    public double[] toPrimitive(ArrayList<Double> in){
        double[] arr = new double[in.size()];

        for(int i = 0; i < in.size(); i++) {
            if (in.get(i) != null) {
                arr[i] = in.get(i);
            }
        }

        return arr;
    }
    
}