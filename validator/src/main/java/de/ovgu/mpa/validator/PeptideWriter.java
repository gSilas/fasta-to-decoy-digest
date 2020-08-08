package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class PeptideWriter {

	List<String[]> aaLines;

	int[] peptideLengthArray;
	int[] peptideLengthRedundantArray;
	HashMap<Integer, Integer> proteinLengthMap;

	public PeptideWriter() {
		aaLines = new ArrayList<String[]>();
		peptideLengthArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];
		peptideLengthRedundantArray = new int[ValidatorConfig.MAXIMUM_PEP_LENGTH];
		proteinLengthMap = new HashMap<Integer, Integer>();

		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++) {
			peptideLengthArray[i] = 0;
			peptideLengthRedundantArray[i] = 0;
		}
	}

	public int readFasta(int batchSize, String batchFiles, File fasta) throws IOException {
		double masswater = Constants.MASS_WATER;

		int batchcount = 1;
		int peptideCount = 0;
		int proteinCount = 0;
		double peptideMass;

		ArrayList<Peptide> databasePeptideList = new ArrayList<Peptide>();
		FASTAFileReader fr = new FASTAFileReader(fasta);
		fr.open();
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(batchFiles + batchcount + ".pep")));
		while (fr.hasNext()) {
			Protein fastaProtein = fr.next();

			if (!proteinLengthMap.containsKey(fastaProtein.sequence.length())) {
				proteinLengthMap.put(fastaProtein.sequence.length(), 1);
			} else {
				int oldVal = proteinLengthMap.get(fastaProtein.sequence.length());
				proteinLengthMap.put(fastaProtein.sequence.length(), oldVal + 1);
			}
			proteinCount++;

			LinkedList<String> peptideList = fastaProtein.getPeptides();
			for (String peptideString : peptideList) {

				peptideCount++;
				if (peptideCount % batchSize == 0) {
					System.out.println(peptideCount + " peptides processed");
					batchcount++;
					// this step is crucial, it sorts the peptides by their mass enabling easy
					// merging later
					databasePeptideList.sort(Peptide.getComparator());
					for (Peptide peptide : databasePeptideList) {
						bw.write(peptide.sequence + ";" + Double.toString(peptide.mass) + "\n");
					}
					databasePeptideList.clear();
					bw.flush();
					bw.close();
					bw = new BufferedWriter(new FileWriter(new File(batchFiles + batchcount + ".pep")));
				}

				peptideMass = 0.0;
				for (char c : peptideString.toCharArray()) {
					AminoAcid aa = AminoAcid.fromOneLetter(c);
					peptideMass += aa.getMass();
				}
				peptideMass -= masswater * (peptideString.length() - 1);
				databasePeptideList.add(new Peptide(peptideString, peptideMass));
			}
		}
		fr.close();

		// handle last incomplete batch
		databasePeptideList.sort(Peptide.getComparator());
		for (Peptide peptide : databasePeptideList) {
			bw.write(peptide.sequence + ";" + Double.toString(peptide.mass) + "\n");
		}
		databasePeptideList.clear();
		bw.flush();
		bw.close();
		System.out.println("Batchfiles, total peptides processed: " + peptideCount);

		return proteinCount;
	}

	public void mergeFiles(int batchSize, String batchFolder, String redundantPepDB) throws IOException {
		System.out.println("Merging Petpide Batch files");
		List<String> result = Files.walk(Paths.get(batchFolder)).map(x -> x.toString()).filter(f -> f.endsWith(".pep"))
				.collect(Collectors.toList());
		BufferedReader[] pepFiles = new BufferedReader[result.size()];
		boolean[] endOfFile = new boolean[result.size()];
		String[] currentPepSequences = new String[result.size()];
		double[] currentPepMasses = new double[result.size()];
		// open readers
		int i = 0;
		for (String fileName : result) {
			pepFiles[i] = new BufferedReader(new FileReader(new File(fileName)));
			String[] split = pepFiles[i].readLine().split(";");
			currentPepSequences[i] = split[0];
			currentPepMasses[i] = Double.valueOf(split[1]);
			i++;
		}
		// open writer, write merge file
		int peptideCount = 0;
		int smallestIndex = -1;
		double smallestMass = -1.0;
		int endOfFileCount = 0;
		BufferedWriter brComplete = new BufferedWriter(new FileWriter(new File(redundantPepDB)));
		while (true) {
			smallestIndex = 0;
			smallestMass = currentPepMasses[0];
			for (int j = 1; j < result.size(); j++) {
				if (!endOfFile[j]) {
					if (currentPepMasses[j] < smallestMass) {
						smallestMass = currentPepMasses[j];
						smallestIndex = j;
					}
				}
			}
			brComplete.write(
					currentPepSequences[smallestIndex] + ";" + Double.toString(currentPepMasses[smallestIndex]) + "\n");
			peptideCount++;
			if (peptideCount % batchSize == 0) {
				System.out.println(peptideCount + " peptides processed");
			}
			String line = pepFiles[smallestIndex].readLine();
			if (line != null) {
				String[] split = line.split(";");
				currentPepSequences[smallestIndex] = split[0];
				currentPepMasses[smallestIndex] = Double.valueOf(split[1]);
			} else {
				endOfFile[smallestIndex] = true;
				endOfFileCount++;
				if (endOfFileCount == result.size()) {
					break;
				}
			}
		}
		brComplete.flush();
		brComplete.close();
		System.out.println("Merging, total peptides processed: " + peptideCount);
	}

	public int removeDuplicates(int batchSize, String redundantPepDB, String nonRedundantPepDB) throws IOException {
		int nonRedundantPeptideCount = 0;
		this.aaLines.clear();
		System.out.println("Creating non redundant petpide database");
		BufferedReader brR = new BufferedReader(new FileReader(new File(redundantPepDB)));
		BufferedWriter bwNR = new BufferedWriter(new FileWriter(new File(nonRedundantPepDB)));
		// KEY = peptide sequence, VALUE = line
		HashMap<String, String> currentPeptides = new HashMap<String, String>();
		HashMap<String, Integer> lengthMap = new HashMap<String, Integer>();
		HashMap<String, Integer> aaMap = new HashMap<String, Integer>();

		String sequence = "";
		double mass = -1.0;
		double currentMass = -1.0;
		String line = brR.readLine();
		while (line != null) {
			String[] split = line.split(";");
			sequence = split[0];
			mass = Double.valueOf(split[1]);
			if (currentMass == -1.0) {
				currentMass = mass;
			} else if (mass != currentMass) {
				currentMass = mass;
				// remove duplicates, write entries
				for (final Map.Entry<String, String> entry : currentPeptides.entrySet()) {
					bwNR.write(entry.getValue() + ";" + lengthMap.get(entry.getKey()) + "\n");
					for (String c : entry.getKey().split("")) {
						if (!aaMap.containsKey(c)) {
							aaMap.put(c, 1);
						} else {
							int oldVal = aaMap.get(c);
							aaMap.put(c, oldVal + 1);
						}
					}

					this.peptideLengthArray[sequence.length()]++;
					this.peptideLengthRedundantArray[sequence.length()]++;

					nonRedundantPeptideCount++;
					if (nonRedundantPeptideCount % batchSize == 0) {
						System.out.println(nonRedundantPeptideCount + " peptides written");
					}
				}
				currentPeptides.clear();
			} else {
				this.peptideLengthRedundantArray[sequence.length()]++;
			}
			if (!currentPeptides.containsKey(sequence)) {
				currentPeptides.put(sequence, line);
				lengthMap.put(sequence, 1);
			} else {
				int oldVal = lengthMap.get(sequence);
				lengthMap.put(sequence, oldVal + 1);
			}
			line = brR.readLine();
		}
		brR.close();
		bwNR.flush();
		bwNR.close();

		aaLines.add(new String[] { "AA", "num" });

		for (final Map.Entry<String, Integer> entry : aaMap.entrySet()) {
			aaLines.add(new String[] { entry.getKey().toString(), entry.getValue().toString() });
		}

		System.out.println("Redundant petpides removed, " + nonRedundantPeptideCount + " peptides written");
		return nonRedundantPeptideCount;
	}

	public void createFiles(String folder, String batchFolder, File fasta) throws IOException {
		// PARAMETERS
		int batchSize = ValidatorConfig.PEPTIDE_BATCH_SIZE;
		String batchFiles = batchFolder + "/Batch_";
		String redundantPepDB = folder + "/Redundant.pep";
		String nonRedundantPepDB = folder + "/NonRedundant.pep";

		System.out.println("Creating Petpide Batch files");
		//
		// STEP ONE
		//
		// fastareader --> create peptides batchfiles
		int proteinCount = readFasta(batchSize, batchFiles, fasta);

		//
		// STEP TWO
		//
		// read batch files --> create merged file
		// sort by mass descending
		mergeFiles(batchSize, batchFolder, redundantPepDB);

		//
		// STEP THREE
		//
		// create duplicate free file
		int peptideCount = removeDuplicates(batchSize, redundantPepDB, nonRedundantPepDB);

		File f = new File(redundantPepDB);
		if (f.delete()) {
			System.out.println(f.getName() + " deleted");
		}

		ArrayList<String[]> peptideLengthLines = new ArrayList<String[]>();
		peptideLengthLines.add(new String[] { "peptide length", "NR occurences", "R occurences" });

		for (int i = 0; i < ValidatorConfig.MAXIMUM_PEP_LENGTH; i++) {
			peptideLengthLines.add(new String[] { Integer.toString(i), Integer.toString(peptideLengthArray[i]), Integer.toString(peptideLengthRedundantArray[i]) });
		}

		ArrayList<String[]> proteinLengthLines = new ArrayList<String[]>();
		proteinLengthLines.add(new String[] { "protein length", "occurences" });
		for (final Map.Entry<Integer, Integer> entry : proteinLengthMap.entrySet()) {
			proteinLengthLines.add(new String[] { entry.getKey().toString(), entry.getValue().toString() });
		}

		final List<String[]> statLines = new ArrayList<>();
		statLines.add(new String[] { "prop", "count" });
		statLines.add(new String[] { "protein count", Integer.toString(proteinCount) });
		statLines.add(new String[] { "peptide count", Integer.toString(peptideCount) });

		final CSVWriter writer = new CSVWriter();
		try {
			writer.createCSV(aaLines, folder + "/aa_num.csv");
			writer.createCSV(statLines, folder + "/stats.csv");
			writer.createCSV(peptideLengthLines, folder + "/peptidelength.csv");
			writer.createCSV(proteinLengthLines, folder + "/proteinlength.csv");
		} catch (final FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		deleteDir(new File(batchFolder));
	}

	void deleteDir(File file) {
		File[] contents = file.listFiles();
		if (contents != null) {
			for (File f : contents) {
				if (!Files.isSymbolicLink(f.toPath())) {
					deleteDir(f);
				}
			}
		}
		file.delete();
	}

}
