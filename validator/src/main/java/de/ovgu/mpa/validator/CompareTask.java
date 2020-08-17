package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class CompareTask implements Callable<CompareResult> {

	ExecutorService pool;

	double[] sqrtLookup;

	final String decoyPeptides; 
	final String targetPeptides; 
	CompareResult result;

	public CompareTask(String decoyPeptides, String targetPeptides) {
		this.decoyPeptides = decoyPeptides;
		this.targetPeptides = targetPeptides;
		this.pool = Executors.newCachedThreadPool();
		sqrtLookup = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH * 4; i++) {
			sqrtLookup[i] = Math.sqrt(i);
		}
		this.result = new CompareResult();
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

		final BufferedReader brD = new BufferedReader(new FileReader(new File(decoyPeptides)));
		final BufferedReader brT = new BufferedReader(new FileReader(new File(targetPeptides)));
		final double tolerance = 0.1;
		final double greatMatchTolerance = 0.5;

		String targetLine = brT.readLine();
		String[] targetSplit = targetLine.split(";");
		String targetSequence = targetSplit[0];
		Double targetMass = Double.valueOf(targetSplit[1]);

		result.lengthTargetArray[targetSequence.length()]++;

		String decoyLine = brD.readLine();
		String[] decoySplit = decoyLine.split(";");
		String decoySequence = decoySplit[0];
		Double decoyMass = Double.valueOf(decoySplit[1]);

		result.lengthDecoyArray[decoySequence.length()]++;

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
				result.lengthDecoyArray[decoySequence.length()]++;

			} else if (targetMass > decoyMass + tolerance) {
				// target larger than decoy
				decoyLine = brD.readLine();
				if (decoyLine == null)
					break;

				decoySplit = decoyLine.split(";");
				decoySequence = decoySplit[0];
				decoyMass = Double.valueOf(decoySplit[1]);
				result.lengthDecoyArray[decoySequence.length()]++;

			} else if (targetMass < decoyMass - tolerance) {
				// target smaller than decoy

				// calculate matches for decoys in window
				for (Future<double[]> decoyIonsFuture : asyncDecoyIonsList) {
					double cosineSimilarity = matchIons(decoyIonsFuture, targetIonsFuture, tolerance);
					result.matchCounter++;

					if (cosineSimilarity >= greatMatchTolerance) {
						result.lengthMatchArray[targetSequence.length()]++;
						result.greatMatchCounter++;
					}

					int bin = (int) (cosineSimilarity * 10);
					result.cosineSimilarityBins[bin]++;
					result.lengthBinMatrix[targetSequence.length()][bin]++;
				}

				// get new target
				targetLine = brT.readLine();

				if (targetLine == null)
					break;

				targetSplit = targetLine.split(";");
				targetSequence = targetSplit[0];
				targetMass = Double.valueOf(targetSplit[1]);
				result.lengthTargetArray[targetSequence.length()]++;

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
					result.matchCounter++;

					if (cosineSimilarity >= greatMatchTolerance) {
						result.lengthMatchArray[targetSequence.length()]++;
						result.greatMatchCounter++;
					}

					int bin = (int) (cosineSimilarity * 10);
					result.cosineSimilarityBins[bin]++;
					result.lengthBinMatrix[targetSequence.length()][bin]++;
				}

				targetLine = brT.readLine();

				if (targetLine == null)
					break;

				targetSplit = targetLine.split(";");
				targetSequence = targetSplit[0];
				targetMass = Double.valueOf(targetSplit[1]);
				result.lengthTargetArray[targetSequence.length()]++;

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

		brT.close();
		brD.close();
	}

	@Override
	public CompareResult call() throws Exception {
		try {
			this.comparePeptides();
		} catch (IOException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
		}
		return result;
	}

}