package de.ovgu.mpa.validator.mascot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Stream;

import de.ovgu.mpa.validator.FragmentationCallable;
import de.ovgu.mpa.validator.ValidatorConfig;
import de.ovgu.mpa.validator.XTandemMatch;

public class MascotCompareCallable implements Callable<MascotCompareResult> {

	private class MatchResult {

		double cosineSimilarity;
		int fragmentMatch;

		public MatchResult(double cosineSimilarity, int fragmentMatch) {
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

	public MatchResult matchIons(Future<double[]> theoreticalIonsFuture, Future<double[][]> experimentalIonsFuture,
			double tolerance) throws InterruptedException, ExecutionException {

		double[] theoreticalIons = theoreticalIonsFuture.get();
		double[] experimentalIons = experimentalIonsFuture.get()[0];

		int indexTheo = 0;
		int indexExp = 0;
		int fragmentMatch = 0;

		while (true) {
			if (experimentalIons[indexExp] > theoreticalIons[indexTheo] + tolerance) {
				indexTheo++;
			} else if (experimentalIons[indexExp] < theoreticalIons[indexTheo] - tolerance) {
				indexExp++;
			} else {
				fragmentMatch++;
				indexExp++;
				indexTheo++;
			}

			if (indexTheo >= theoreticalIons.length || indexExp >= experimentalIons.length) {
				break;
			}
		}

		double cosineSimilarity = fragmentMatch
				/ (Math.sqrt(experimentalIons.length) * sqrtLookup[theoreticalIons.length]);
		return new MatchResult(cosineSimilarity, fragmentMatch);
	}

	public Future<double[][]> parseIons(String mascotIonsString) {
		return pool.submit(() -> {
			double[][] mascotResult = new double[2][mascotIonsString.split(":").length - 1];
			double[] mascotIons = new double[mascotIonsString.split(":").length - 1];
			Double[] mascotIntensity = new Double[mascotIonsString.split(":").length - 1];

			int i = 0;
			for (String ionPair : mascotIonsString.split(",")) {
				mascotIons[i] = Double.parseDouble(ionPair.split(":")[0]);
				mascotIntensity[i] = Double.parseDouble(ionPair.split(":")[1]);
				i++;
			}

			final List<Double> mascotIntensityCopy = Arrays.asList(mascotIntensity);
			ArrayList<Double> sortedList = new ArrayList(mascotIntensityCopy);
			Collections.sort(sortedList, Comparator.comparing(s -> mascotIons[mascotIntensityCopy.indexOf(s)]));
			Arrays.sort(mascotIons);

			mascotResult[0] = mascotIons;
			mascotResult[1] = new double[mascotIonsString.split(":").length - 1];

			int j = 0;
			for (Double d : sortedList) {
				mascotResult[1][j] = (double) d;
				j++;
			}

			return mascotResult;
		});
	}

	public void comparePeptides() throws IOException, InterruptedException, ExecutionException {

		final BufferedReader brM = new BufferedReader(new FileReader(new File(mascotPeptides)));
		final double tolerance = 0.005;

		String mascotLine = brM.readLine();
		String[] mascotSplit = mascotLine.split(";");
		String mascotSequence = mascotSplit[0];
		// Double mascotMass = Double.valueOf(mascotSplit[1]);
		Double mascotScore = Double.valueOf(mascotSplit[3]);
		String mascotIonsString = mascotSplit[4];

		Future<double[][]> mascotFuture = parseIons(mascotIonsString);
		Future<double[]> mascotIonsSeqFuture = pool.submit(new FragmentationCallable(mascotSequence, true));

		while (true) {

			MatchResult match = matchIons(mascotIonsSeqFuture, mascotFuture, tolerance);
			XTandemMatch xTandemMatch = new XTandemMatch(mascotSequence, mascotFuture, pool);
			double xTandemScore = xTandemMatch.calculate();

			result.cosineSimilarityMap.put(mascotSequence, match.cosineSimilarity);
			result.mascotScoreMap.put(mascotSequence, mascotScore);
			result.fragmentMatchMap.put(mascotSequence, match.fragmentMatch);
			result.xTandemMap.put(mascotSequence, xTandemScore);

			mascotLine = brM.readLine();

			if (mascotLine == null) {
				break;
			}

			mascotSplit = mascotLine.split(";");
			mascotSequence = mascotSplit[0];
			// mascotMass = Double.valueOf(mascotSplit[1]);
			mascotScore = Double.valueOf(mascotSplit[3]);
			mascotIonsString = mascotSplit[4];

			mascotFuture = parseIons(mascotIonsString);
			mascotIonsSeqFuture = pool.submit(new FragmentationCallable(mascotSequence, true));

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