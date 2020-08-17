package de.ovgu.mpa.validator;

import static java.nio.file.StandardCopyOption.*;
import org.apache.commons.cli.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class App {

	public static void main(String[] args) {

		Options options = new Options();

		Option fasta = new Option("f", "read_fasta", true, "fasta file path");
		fasta.setRequired(false);
		options.addOption(fasta);

		Option compare = new Option("c", "compare_dbs", true, "paths to db1 abd db2");
		compare.setRequired(false);
		options.addOption(compare);

		Option excludeX = new Option("ex", "exclude_x", false,
				"exclude peptides with unkown amino acids; requires read_fasta");
				excludeX.setRequired(false);
		options.addOption(excludeX);

		Option timer = new Option("t", "timer", false, "time execution");
		timer.setRequired(false);
		options.addOption(timer);

		Option threads = new Option("p", "threads", true, "number of thread");
		threads.setRequired(false);
		options.addOption(threads);

		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd = null;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("Fasta DB Compare", options);
			System.exit(1);
		}

		if (cmd.hasOption("exclude_x")) {
			System.out.println("Excluding X");
			ValidatorConfig.excludeX = true;
		}

		int numThreads = 1;
		if (cmd.hasOption("threads")) {
			try {
				System.out.println(cmd.getOptionValue("threads"));
				numThreads = Integer.parseInt(cmd.getOptionValue("threads"));
			} catch (NumberFormatException e) {
				System.err.println("Could not parse threads option, number of threads set to one.");
				numThreads = 1;
			}
		}

		long startTime = System.currentTimeMillis();

		if (cmd.hasOption("read_fasta")) {
			File batchDir = new File("tmp_batches");
			if (!batchDir.exists())
				batchDir.mkdir();
			File fastaFolder = new File("fasta");
			if (!fastaFolder.exists())
				fastaFolder.mkdir();

			// Create Target and Decoy Folder
			String fileName = cmd.getOptionValue("read_fasta")
					.split("/")[cmd.getOptionValue("read_fasta").split("/").length - 1].split("\\.")[0];
			File targetFolder = new File(fastaFolder.getPath() + "/" + fileName);
			if (!targetFolder.exists())
				targetFolder.mkdir();
			File decoyFolder = new File(fastaFolder.getPath() + "/" + fileName + "_decoy");
			if (!decoyFolder.exists())
				decoyFolder.mkdir();

			processFasta(targetFolder.toString(), decoyFolder.toString(), cmd.getOptionValue("read_fasta"),
					batchDir.toString(), fastaFolder);
		} else if (cmd.hasOption("compare_dbs")) {
			File resultsFolder = new File("results");
			if (!resultsFolder.exists())
				resultsFolder.mkdir();

			String targetFolder = cmd.getOptionValue("compare_dbs").split(" ")[0];
			String decoyFolder = cmd.getOptionValue("compare_dbs").split(" ")[1];
			compareDB(targetFolder.toString(), decoyFolder.toString(), resultsFolder.toString(), numThreads);

		} else {
			formatter.printHelp("Fasta DB Compare", options);
			System.exit(1);
		}

		if (cmd.hasOption("timer")) {
			System.out.println("execution time: " + ((double) System.currentTimeMillis() - (double) startTime) / 1000.0
					+ " seconds");
		}
	}

	public static void processFasta(String targetFolder, String decoyFolder, String fastaPath, String batchFolder,
			File fastaFolder) {
		// Copy Fasta
		Path source = Paths.get(fastaPath);
		Path target = Paths.get(targetFolder + "/" + fastaPath.split("/")[fastaPath.split("/").length - 1]);
		try {
			Files.copy(source, target, REPLACE_EXISTING);
		} catch (IOException e2) {
			e2.printStackTrace();
		}
		// process fasta
		PeptideWriter writer = new PeptideWriter();
		try {
			writer.createFiles(targetFolder, batchFolder, target.toFile());
		} catch (IOException e) {
			e.printStackTrace();
		}

		// create decoy
		DecoyGenerator.generateReverseDecoyDatabase(target,
				Paths.get(decoyFolder + "/decoy_" + fastaPath.split("/")[fastaPath.split("/").length - 1]));

		writer = new PeptideWriter();
		File batchDir = new File("tmp_batches");
		if (!batchDir.exists())
			batchDir.mkdir();

		try {
			writer.createFiles(decoyFolder, batchFolder, Paths
					.get(decoyFolder + "/decoy_" + fastaPath.split("/")[fastaPath.split("/").length - 1]).toFile());
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void compareDB(String targetFolder, String decoyFolder, String resultsFolder, int numThreads) {
		System.out.println("compare");
		// prepare database splits and thread pool
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads);
		String threadFolder = resultsFolder + "/threads/";
		if (!(new File(threadFolder)).exists())
			(new File(threadFolder)).mkdir();

		HashSet<String> threadNames = new HashSet<String>();
		HashMap<String, String> targetFolders = new HashMap<>();

		// count entries, TODO: replace with readout of peptide count from stats
		long targetCount = 0;
		try {
			BufferedReader brCounting = new BufferedReader(
					new FileReader(new File(targetFolder)));
			String line = brCounting.readLine();
			while (line != null) {
				targetCount++;
				line = brCounting.readLine();
			}
			brCounting.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(targetCount);
		// create subfiles based on number of threads
		try {
			long pepsPerThread = targetCount / numThreads;
			System.out.println(pepsPerThread);
			BufferedReader brDividing = new BufferedReader(
					new FileReader(new File(targetFolder)));
			for (int i = 0; i < numThreads; i++) {
				String threadName = "Thread_" + i;
				threadNames.add(threadName);
				BufferedWriter bwSubfile = new BufferedWriter(
						new FileWriter(new File(threadFolder + "/Target_" + threadName)));
				String line = brDividing.readLine();
				long count = 0;
				while (true) {
					// if also handles last entry
					if (count > pepsPerThread || line == null) {
						if (line != null) {
							bwSubfile.write(line + "\n");
						}
						break;
					} else {
						bwSubfile.write(line + "\n");
						line = brDividing.readLine();
						count++;
					}
				}
				bwSubfile.close();
				targetFolders.put(threadName, threadFolder + "/Target_" + threadName);
				if (!(new File(threadFolder + "/" + threadName)).exists())
					(new File(threadFolder + "/" + threadName)).mkdir();
			}
			brDividing.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		LinkedList<Future<CompareResult>> resultList = new LinkedList<Future<CompareResult>>();

		for (String thread : threadNames) {
			CompareTask task = new CompareTask(decoyFolder, targetFolders.get(thread));
			resultList.add(threadPool.submit(task));
		}
		threadPool.shutdown();

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

		long matchCounter = 0;
		long greatMatchCounter = 0;

		for (Future<CompareResult> resultFuture : resultList) {
			CompareResult result;
			try {
				result = resultFuture.get();

				for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++) {
					lengthTargetArray[i] += result.lengthTargetArray[i];
					lengthDecoyArray[i] += result.lengthDecoyArray[i];
					lengthMatchArray[i] += result.lengthMatchArray[i];
	
					for (int j = 0; j < 11; j++){
						lengthBinMatrix[i][j] += result.lengthBinMatrix[i][j];
					}
				}
	
				for (int j = 0; j < 11; j++){
					cosineSimilarityBins[j] += result.cosineSimilarityBins[j];
				}
	
				matchCounter += result.matchCounter;
				greatMatchCounter += result.greatMatchCounter;

			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		System.out.println("Recognized " + matchCounter + " matches!");
		System.out.println("Recognized " + greatMatchCounter + " great matches!");

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

		final List<String[]> cosineSimilarityBinsLines = new ArrayList<>();
		cosineSimilarityBinsLines.add(new String[] { "bin", "occurences" });
		for (int i = 0; i < cosineSimilarityBins.length; i++) {
			cosineSimilarityBinsLines
			.add(new String[] { Double.toString((double) i / 10.0), Long.toString(cosineSimilarityBins[i]) });
		}

		// combine results
		final CSVWriter writer = new CSVWriter();
		try {
			writer.createCSV(lengthMatchLines, resultsFolder + "/" + "lengthMatch.csv");
			writer.createCSV(cosineSimilarityBinsLines, resultsFolder + "/" + "cosineSimilarityBins.csv");
			writer.createCSV(lengthCollectionLines, resultsFolder + "/" + "lengthCollection.csv");
		} catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// TODO
	}
}
