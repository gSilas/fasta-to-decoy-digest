package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Statistics {

    private final double[] sqrtLookup;

    private double[] targetIons;
    private double[] decoyIons;

    private HashMap<String, double[]> decoyIonMap;

    public Statistics() {
        sqrtLookup = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
        for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH * 4; i++) {
            sqrtLookup[i] = Math.sqrt(i);
        }

        targetIons = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
        decoyIons = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
        decoyIonMap = new HashMap<String, double[]>();
    }

    private double matchPeptides(double tolerance, String targetSequence, String decoySequence,
            HashMap<Integer, Integer> fragmentMatchMap) {
        // MATCH
        FragmentationCalculator.getFragmentIons(targetSequence, targetIons);
        
        if (!decoyIonMap.containsKey(decoySequence)) {
            FragmentationCalculator.getFragmentIons(decoySequence, decoyIons);
            decoyIonMap.put(decoySequence, decoyIons.clone());
        } else {
            decoyIons = decoyIonMap.get(decoySequence);
        }

        int indexDecoy = 0;
        int indexTarget = 0;
        int fragmentMatch = 0;

        double decoyIon = decoyIons[0];
        double targetIon = targetIons[0];

        while (true) {
            if (targetIon <= decoyIon + tolerance && targetIon >= decoyIon - tolerance) {
                fragmentMatch++;
                indexTarget++;
                indexDecoy++;

                if (indexDecoy >= decoySequence.length() || indexTarget >= targetSequence.length()) {
                    break;
                }

                decoyIon = decoyIons[indexDecoy];
                targetIon = targetIons[indexTarget];

            } else if (targetIon > decoyIon + tolerance) {
                indexDecoy++;

                if (indexDecoy >= decoySequence.length()) {
                    break;
                }

                decoyIon = decoyIons[indexDecoy];
            } else if (targetIon < decoyIon - tolerance) {
                indexTarget++;

                if (indexTarget >= targetSequence.length()) {
                    break;
                }

                targetIon = targetIons[indexTarget];
            }
        }

        double cosineSimilarity = fragmentMatch
                / (sqrtLookup[targetSequence.length()] * sqrtLookup[decoySequence.length()]);

        if (fragmentMatchMap.containsKey(fragmentMatch)) {
            fragmentMatchMap.put(fragmentMatch, fragmentMatchMap.get(fragmentMatch) + 1);
        } else {
            fragmentMatchMap.put(fragmentMatch, 1);
        }

        return cosineSimilarity;
    }

    public void comparePeptides(final String decoyPeptides, final String targetPeptides, String folder)
            throws IOException {
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
        int[][] lengthBinMatrix = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH][11];

        // <length, great matches>
        // HashMap<Integer, Integer> lengthMap = new HashMap<>();
        int[] lengthMatchArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];

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

        // <FragmentMatches, Occurences>
        HashMap<Integer, Integer> fragmentMatchMap = new HashMap<>();

        int[] cosineSimilarityBins = new int[11];

        LinkedList<String> decoySequenceWindow = new LinkedList<String>();
        LinkedList<Double> decoyMassWindow = new LinkedList<Double>();

        while (true) {

            if (targetMass <= decoyMass + tolerance && targetMass >= decoyMass - tolerance) {
                // decoy matches target

                decoySequenceWindow.add(decoySequence);
                decoyMassWindow.add(decoyMass);

                decoyLine = brD.readLine();

                if (decoyLine == null)
                    break;

                decoySplit = decoyLine.split(";");
                // splitPeptideLine(decoyLine, ';', decoySplit);
                decoySequence = decoySplit[0];
                decoyMass = Double.valueOf(decoySplit[1]);
                lengthDecoyArray[decoySequence.length()]++;

            } else if (targetMass > decoyMass + tolerance) {
                // target larger than decoy
                decoyLine = brD.readLine();

                if (decoyLine == null)
                    break;

                decoySplit = decoyLine.split(";");
                // splitPeptideLine(decoyLine, ';', decoySplit);
                decoySequence = decoySplit[0];
                decoyMass = Double.valueOf(decoySplit[1]);
                lengthDecoyArray[decoySequence.length()]++;

            } else if (targetMass < decoyMass - tolerance) {
                // target smaller than decoy

                // calculate matches for decoys in window
                for (int i = 0; i < decoySequenceWindow.size(); i++) {
                    String sequence = decoySequenceWindow.get(i);
                    double cosineSimilarity = matchPeptides(tolerance, targetSequence, sequence, fragmentMatchMap);

                    matchCounter++;

                    if (cosineSimilarity >= greatMatchTolerance) {
                        lengthMatchArray[targetSequence.length()]++;
                        greatMatchCounter++;
                    }

                    int bin = (int) (cosineSimilarity * 10);
                    cosineSimilarityBins[bin]++;
                    lengthBinMatrix[targetSequence.length()][bin]++;
                }

                // get new target
                targetLine = brT.readLine();

                if (targetLine == null)
                    break;

                targetSplit = targetLine.split(";");
                // splitPeptideLine(targetLine, ';', targetSplit);
                targetSequence = targetSplit[0];
                targetMass = Double.valueOf(targetSplit[1]);
                lengthTargetArray[targetSequence.length()]++;

                // check if decoys in window don't match with target and remove unmatching
                // decoys

                while (!decoySequenceWindow.isEmpty()) {
                    double mass = decoyMassWindow.peek();
                    if (targetMass > mass + tolerance) {
                        // target larger than decoy
                        decoyMassWindow.pop();
                        String sequence = decoySequenceWindow.pop();
                        decoyIonMap.remove(sequence);
                    } else {
                        break;
                    }
                }
            }
        }

        if (decoyLine == null) {
            while (!decoyMassWindow.isEmpty()) {
                for (int i = 0; i < decoySequenceWindow.size(); i++) {
                    String sequence = decoySequenceWindow.get(i);
                    double cosineSimilarity = matchPeptides(tolerance, targetSequence, sequence, fragmentMatchMap);

                    matchCounter++;

                    if (cosineSimilarity >= greatMatchTolerance) {
                        lengthMatchArray[targetSequence.length()]++;
                        greatMatchCounter++;
                    }

                    int bin = (int) (cosineSimilarity * 10);
                    cosineSimilarityBins[bin]++;
                    lengthBinMatrix[targetSequence.length()][bin]++;

                }

                targetLine = brT.readLine();

                if (targetLine == null)
                    break;

                targetSplit = targetLine.split(";");
                // splitPeptideLine(targetLine, ';', targetSplit);
                targetSequence = targetSplit[0];
                targetMass = Double.valueOf(targetSplit[1]);
                lengthTargetArray[targetSequence.length()]++;

                // check if decoys in window don't match with target and remove unmatching
                // decoys

                while (!decoySequenceWindow.isEmpty()) {
                    double mass = decoyMassWindow.peek();
                    if (targetMass > mass + tolerance) {
                        // target larger than decoy
                        decoyMassWindow.pop();
                        String sequence = decoySequenceWindow.pop();
                        decoyIonMap.remove(sequence);
                    } else {
                        break;
                    }
                }
            }
        }

        final List<String[]> lengthMatchLines = new ArrayList<>();
        final List<String[]> lengthCollectionLines = new ArrayList<>();

        lengthMatchLines.add(new String[] { "length", "matches" });
        lengthCollectionLines.add(new String[] { "length", "target", "decoy", "matches", "bins" });

        for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++) {
            lengthMatchLines.add(new String[] { Integer.toString(i), Integer.toString(lengthMatchArray[i]) });

            String bString = "";
            for (int j = 0; j < 10; j++) {
                bString += "," + lengthBinMatrix[i][j];
            }
            lengthCollectionLines.add(new String[] { Integer.toString(i), Integer.toString(lengthTargetArray[i]),
                    Integer.toString(lengthDecoyArray[i]), Integer.toString(lengthMatchArray[i]), bString });
        }
        final List<String[]> fragmentMatchLines = new ArrayList<>();
        fragmentMatchLines.add(new String[] { "matches", "relative occurences" });

        for (final Map.Entry<Integer, Integer> entry : fragmentMatchMap.entrySet()) {
            fragmentMatchLines.add(new String[] { entry.getKey().toString(),
                    Double.toString(100 * ((double) entry.getValue() / (double) matchCounter)) });
        }

        final List<String[]> cosineSimilarityBinsLines = new ArrayList<>();
        cosineSimilarityBinsLines.add(new String[] { "bin", "occurences" });
        for (int i = 0; i < cosineSimilarityBins.length; i++) {
            cosineSimilarityBinsLines
                    .add(new String[] { Double.toString((double) i / 10.0), Integer.toString(cosineSimilarityBins[i]) });
        }

        final CSVWriter writer = new CSVWriter();
        try {
            writer.createCSV(fragmentMatchLines, folder + "/" + "fragmentMatch.csv");
            writer.createCSV(lengthMatchLines, folder + "/" + "lengthMatch.csv");
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

}