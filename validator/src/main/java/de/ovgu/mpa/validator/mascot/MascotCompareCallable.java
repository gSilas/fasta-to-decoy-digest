package de.ovgu.mpa.validator.mascot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;

import de.ovgu.mpa.validator.FragmentationCallable;
import de.ovgu.mpa.validator.ValidatorConfig;

public class MascotCompareCallable implements Callable<MascotCompareResult> {

	private class MatchResult {

		double cosineSimilarity;
		int fragmentMatch;

		public MatchResult(double cosineSimilarity, int fragmentMatch){
			this.cosineSimilarity = cosineSimilarity;
			this.fragmentMatch = fragmentMatch;
		}
	}

	ExecutorService pool;

	double[] sqrtLookup;

	final String mascotPeptides; 
	MascotCompareResult result;
	AtomicLong done;

	public MascotCompareCallable(String mascotPeptides, AtomicLong done) {
		this.mascotPeptides = mascotPeptides;
		this.pool = Executors.newCachedThreadPool();
		sqrtLookup = new double[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];
		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH * 4; i++) {
			sqrtLookup[i] = Math.sqrt(i);
		}
		this.result = new MascotCompareResult();
		this.done = done;
	}

	public MatchResult matchIons(Future<double[]> theoreticalIonsFuture, Future<double[]> experimentalIonsFuture, double tolerance)
			throws InterruptedException, ExecutionException {

		double[] decoyIons = theoreticalIonsFuture.get();
		double[] targetIons = experimentalIonsFuture.get();

		//System.out.println(Arrays.toString(decoyIons));
		//System.out.println(Arrays.toString(targetIons));

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

				// System.out.println("Match");

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

		double cosineSimilarity = fragmentMatch / (sqrtLookup[decoyIons.length] * sqrtLookup[decoyIons.length]);
		return new MatchResult(cosineSimilarity, fragmentMatch);
	}

	public Future<double[]> parseIons(String mascotIonsString) {
		return pool.submit(() -> {
			double[] mascotIons = new double[mascotIonsString.split(":").length];
			// System.out.println(mascotIonsString);
			int i = 0;
			for (String ionPair : mascotIonsString.split(",")) {
				if (Double.parseDouble(ionPair.split(":")[1]) > 25.0) {
					mascotIons[i] = Double.parseDouble(ionPair.split(":")[0]);
					i++;
				}
			}
			Arrays.sort(mascotIons);
			return mascotIons;
		});
	}

	public void comparePeptides()
			throws IOException, InterruptedException, ExecutionException {

		final BufferedReader brM = new BufferedReader(new FileReader(new File(mascotPeptides)));
		final double tolerance = 0.1;
		
		String mascotLine = brM.readLine();
		String[] mascotSplit = mascotLine.split(";");
		String mascotSequence = mascotSplit[0];
		// Double mascotMass = Double.valueOf(mascotSplit[1]);
		Double mascotScore = Double.valueOf(mascotSplit[3]);
		String mascotIonsString = mascotSplit[4];

		Future<double[]> mascotIonsFuture = parseIons(mascotIonsString);

		this.done.incrementAndGet();

		Future<double[]> mascotIonsSeqFuture = pool.submit(new FragmentationCallable(mascotSequence));

		// double bestMatchScore = Double.NEGATIVE_INFINITY;

		while (true) {

			MatchResult match = matchIons(mascotIonsSeqFuture, mascotIonsFuture, tolerance);

			result.cosineSimilarityMap.put(mascotSequence, match.cosineSimilarity);
			result.mascotScoreMap.put(mascotSequence, mascotScore);
			result.fragmentMatchMap.put(mascotSequence, match.fragmentMatch);

			mascotLine = brM.readLine();

			if (mascotLine == null) {
				break;
			}

			mascotSplit = mascotLine.split(";");
			mascotSequence = mascotSplit[0];
			// mascotMass = Double.valueOf(mascotSplit[1]);
			mascotScore = Double.valueOf(mascotSplit[3]);
			mascotIonsString = mascotSplit[4];

			mascotIonsFuture = parseIons(mascotIonsString);

		}

		pool.shutdown();

		brM.close();
	}

	@Override
	public MascotCompareResult call() throws Exception {
		try {
			this.comparePeptides();
		} catch (IOException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
		}
		return result;
	}

}