package de.ovgu.mpa.validator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

public class FASTAFileReader implements Iterator<FastaProtein> {

	private final File fastaFile;
	
	private BufferedReader fastaReader;
	private StringBuilder sequencebuilder;
	private FastaProtein currentProtein;
	private String nextDescription;

	public FASTAFileReader(File fastafile) throws IOException {
		this.fastaFile = fastafile;
	}
	
	public void open() throws FileNotFoundException {
		this.nextDescription = null;
		this.currentProtein = new FastaProtein();
		this.sequencebuilder = new StringBuilder();
		this.fastaReader = new BufferedReader(new FileReader(fastaFile));
	}
	
	public void close() throws IOException {
		this.fastaReader.close();
	}

	@Override
	public FastaProtein next() {
		this.currentProtein.setSequence(this.sequencebuilder.toString().replaceAll("\n", ""));
		return this.currentProtein;
	}

	@Override
	public boolean hasNext() {
		try {
			this.currentProtein = new FastaProtein();
			this.sequencebuilder.setLength(0);
			String line = this.fastaReader.readLine();
			while (line != null) {
				if (line.startsWith(">") && this.nextDescription != null) {
					// next protein
					this.currentProtein.setDescription(nextDescription.substring(1));
					this.nextDescription = line;
					return true;
				} else if (line.startsWith(">") && this.nextDescription == null) {
					// first protein
					this.nextDescription = line;
				} else {
					// sequence
					this.sequencebuilder.append(line);
				}
				line = this.fastaReader.readLine();
				if (line == null) {
					// last protein
					this.currentProtein.setDescription(nextDescription.substring(1));
					return true;
				}
			}
		} catch (IOException e)  {
			e.printStackTrace();
		}
		return false;
	}

}
