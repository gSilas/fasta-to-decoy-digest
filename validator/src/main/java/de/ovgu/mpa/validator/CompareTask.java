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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class CompareTask implements Runnable {

	ExecutorService pool;

	double[] sqrtLookup;

	final String decoyPeptides; 
	final String targetPeptides; 
	String folder;

	public CompareTask(String decoyPeptides, String targetPeptides, String folder) {
		this.decoyPeptides = decoyPeptides;
		this.targetPeptides = targetPeptides;
		this.folder = folder;
		this.pool = Executors.newCachedThreadPool();
		sqrtLookup = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH * 4; i++) {
			sqrtLookup[i] = Math.sqrt(i);
		}

	}

	@Override
	public void run() {
		try {
			this.comparePeptides();
		} catch (IOException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
		}
	}

	public Double matchIons(Future<double[]> decoyIonsFuture, Future<double[]> targetIonsFuture, double tolerance)
			throws InterruptedException, ExecutionException {

		double[] decoyIons = decoyIonsFuture.get();
		double[] targetIons = targetIonsFuture.get();

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

				if (indexDecoy >= decoyIons.length || indexTarget >= targetIons.length) {
					break;
				}

				decoyIon = decoyIons[indexDecoy];
				targetIon = targetIons[indexTarget];

			} else if (targetIon > decoyIon + tolerance) {
				indexDecoy++;

				if (indexDecoy >= decoyIons.length) {
					break;
				}

				decoyIon = decoyIons[indexDecoy];
			} else if (targetIon < decoyIon - tolerance) {
				indexTarget++;

				if (indexTarget >= targetIons.length) {
					break;
				}

				targetIon = targetIons[indexTarget];
			}

		}

		double cosineSimilarity = fragmentMatch / (sqrtLookup[targetIons.length] * sqrtLookup[decoyIons.length]);

		return cosineSimilarity;
	}

	public void comparePeptides()
			throws IOException, InterruptedException, ExecutionException {

		final BufferedReader brD = new BufferedReader(new FileReader(new File(decoyPeptides + "/NonRedundant.pep")));
		final BufferedReader brT = new BufferedReader(new FileReader(new File(targetPeptides)));
		final double tolerance = 0.1;
		final double greatMatchTolerance = 0.5;
		long matchCounter = 0;
		long greatMatchCounter = 0;

		// <length, target proteins>
		long[] lengthTargetArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];

		// <length, decoy proteins>
		long[] lengthDecoyArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];

		// <length, bins>
		long[][] lengthBinMatrix = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH][11];

		// <length, great matches>
		// HashMap<Integer, Integer> lengthMap = new HashMap<>();
		long[] lengthMatchArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];

		long[] cosineSimilarityBins = new long[11];

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

		LinkedList<String> decoySequenceWindow = new LinkedList<String>();
		LinkedList<Double> decoyMassWindow = new LinkedList<Double>();

		LinkedList<Future<double[]>> asyncDecoyIonsList = new LinkedList<Future<double[]>>();
		Future<double[]> targetIonsFuture = pool.submit(new FragmentationCallable(targetSequence));

		while (true) {

			if (targetMass <= decoyMass + tolerance && targetMass >= decoyMass - tolerance) {
				// decoy matches target

				decoySequenceWindow.add(decoySequence);
				decoyMassWindow.add(decoyMass);
				Future<double[]> decoyIonsFuture = pool.submit(new FragmentationCallable(decoySequence));
				asyncDecoyIonsList.add(decoyIonsFuture);

				decoyLine = brD.readLine();

				if (decoyLine == null)
					break;

				decoySplit = decoyLine.split(";");
				decoySequence = decoySplit[0];
				decoyMass = Double.valueOf(decoySplit[1]);
				lengthDecoyArray[decoySequence.length()]++;

			} else if (targetMass > decoyMass + tolerance) {
				// target larger than decoy
				decoyLine = brD.readLine();
				if (decoyLine == null)
					break;

				decoySplit = decoyLine.split(";");
				decoySequence = decoySplit[0];
				decoyMass = Double.valueOf(decoySplit[1]);
				lengthDecoyArray[decoySequence.length()]++;

			} else if (targetMass < decoyMass - tolerance) {
				// target smaller than decoy

				// calculate matches for decoys in window
				for (Future<double[]> decoyIonsFuture : asyncDecoyIonsList) {
					double cosineSimilarity = matchIons(decoyIonsFuture, targetIonsFuture, tolerance);
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
				targetSequence = targetSplit[0];
				targetMass = Double.valueOf(targetSplit[1]);
				lengthTargetArray[targetSequence.length()]++;

				targetIonsFuture = pool.submit(new FragmentationCallable(targetSequence));

				// check if decoys in window don't match with target and remove unmatching
				// decoys

				while (!decoySequenceWindow.isEmpty()) {
					double mass = decoyMassWindow.peek();
					if (targetMass > mass + tolerance) {
						// target larger than decoy
						decoyMassWindow.pop();
						decoySequenceWindow.pop();
						asyncDecoyIonsList.pop();
					} else {
						break;
					}
				}

			}
		}

		if (decoyLine == null) {
			while (!decoyMassWindow.isEmpty()) {

				for (Future<double[]> decoyIonsFuture : asyncDecoyIonsList) {
					double cosineSimilarity = matchIons(decoyIonsFuture, targetIonsFuture, tolerance);
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
				targetSequence = targetSplit[0];
				targetMass = Double.valueOf(targetSplit[1]);
				lengthTargetArray[targetSequence.length()]++;

				targetIonsFuture = pool.submit(new FragmentationCallable(targetSequence));

				// check if decoys in window don't match with target and remove unmatching
				// decoys

				while (!decoySequenceWindow.isEmpty()) {
					double mass = decoyMassWindow.peek();
					if (targetMass > mass + tolerance) {
						// target larger than decoy
						decoyMassWindow.pop();
						decoySequenceWindow.pop();
						asyncDecoyIonsList.pop();
					} else {
						break;
					}
				}
			}
		}

		System.out.println("Finished reading target & decoy lists");

		pool.shutdown();

		final List<String[]> lengthMatchLines = new ArrayList<>();
		final List<String[]> lengthCollectionLines = new ArrayList<>();

		lengthMatchLines.add(new String[] { "length", "matches" });
		lengthCollectionLines.add(new String[] { "length", "target", "decoy", "matches", "bins" });

		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++) {
			lengthMatchLines.add(new String[] { Integer.toString(i), Long.toString(lengthMatchArray[i]) });

			String bString = "";
			for (int j = 0; j < 10; j++) {
				bString += "," + lengthBinMatrix[i][j];
			}
			lengthCollectionLines.add(new String[] { Integer.toString(i), Long.toString(lengthTargetArray[i]),
					Long.toString(lengthDecoyArray[i]), Long.toString(lengthMatchArray[i]), bString });
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
			.add(new String[] { Double.toString((double) i / 10.0), Long.toString(cosineSimilarityBins[i]) });
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