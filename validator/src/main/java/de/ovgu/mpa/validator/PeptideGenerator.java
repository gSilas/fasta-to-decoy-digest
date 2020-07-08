package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

public class PeptideGenerator {
	
	public PeptideGenerator() {
	}
	
	public void createFiles(String folder, String batchFolder, File fasta) throws IOException {
		// PARAMETERS
		int batchSize = ValidatorConfig.PEPTIDE_BATCH_SIZE;
		double masswater = Constants.MASS_WATER;
		String batchFiles = batchFolder + "/Batch_";
		String redundantPepDB = folder + "/Redundant.pep";
		String nonRedundantPepDB = folder + "/NonRedundant.pep";
		
		System.out.println("Creating Petpide Batch files");
		//
		// STEP ONE
		//
		// fastareader --> create peptides batchfiles
		int batchcount = 1;
		int peptideCount = 0;
		double peptideMass;
		ArrayList<Peptide> databasePeptideList = new ArrayList<Peptide>(); 
		FASTAFileReader fr = new FASTAFileReader(fasta);
		fr.open();
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(batchFiles + batchcount + ".pep")));
		while (fr.hasNext()) {
			Protein fastaProtein = fr.next();
			LinkedList<String> peptideList = fastaProtein.getPeptides();
			for (String peptideString : peptideList) {
				peptideCount++;
				if (peptideCount % batchSize == 0) {
					System.out.println(peptideCount + " peptides processed");
					batchcount++;
					// this step is crucial, it sorts the peptides by their mass enabling easy merging later
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
		//
		// STEP TWO
		//
		// read batch files --> create merged file
		// sort by mass descending
		System.out.println("Merging Petpide Batch files");
        List<String> result = Files.walk(Paths.get(batchFolder)).map(x -> x.toString())
				.filter(f -> f.endsWith(".pep")).collect(Collectors.toList());
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
		peptideCount = 0;
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
			brComplete.write(currentPepSequences[smallestIndex] + ";" + Double.toString(currentPepMasses[smallestIndex]) + "\n");
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
		
		//
		// STEP THREE
		//
		// create duplicate free file
		int nonRedundantPeptideCount = 0;
		System.out.println("Creating non redundant petpide database");
		BufferedReader brR = new BufferedReader(new FileReader(new File(redundantPepDB)));  
		BufferedWriter bwNR = new BufferedWriter(new FileWriter(new File(nonRedundantPepDB)));
		// KEY = peptide sequence, VALUE = line 
		HashMap<String, String> currentPeptides = new HashMap<String, String>();
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
				for (String value : currentPeptides.values()) {
					bwNR.write(value + "\n");
					nonRedundantPeptideCount++;
					if (nonRedundantPeptideCount % batchSize == 0) {
						System.out.println(nonRedundantPeptideCount + " peptides written");
					}
				}
				currentPeptides.clear();
			}
			if (!currentPeptides.containsKey(sequence)) {
				currentPeptides.put(sequence, line);
			}
			line = brR.readLine();
		}
		brR.close();
		bwNR.flush();
		bwNR.close();
		System.out.println("Redundant petpides removed, " + nonRedundantPeptideCount + " of " + peptideCount + " written");
	}

}

