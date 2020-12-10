package de.ovgu.mpa.validator;

import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class XTandemMatch {

	private class Score {
		final double score;
		final int matchedIons;
		double matchedIonIntensity;
		double matchedIonToTotal;
		int matchedIonNumLongestConsecutiveMatch;
		double[] matchErrorListIon;
		double sumOfMatchedIntensities;

		public Score(double score, int matchedIons) {
			this.score = score;
			this.matchedIons = matchedIons;
		}
	}

	private final double PROTON = 1.007276;

	final double[] pepMZ;
	final double[] pepINT;

	final String sequence;

	ExecutorService pool;

	public XTandemMatch(String sequence, Future<double[][]> pepResult, ExecutorService pool) throws InterruptedException, ExecutionException {
		this.pepMZ = pepResult.get()[0];
		this.pepINT = pepResult.get()[1];

		this.sequence = sequence;
		this.pool = pool;
	}

	public Score score(double[] pepMZDouble, double[] pepINT, double[] specMZDouble, double[] specINT) {

		long[] specMZ = new long[specMZDouble.length];
		for (int i = 0; i < specMZDouble.length; i++){
			specMZ[i] = Math.round(specMZDouble[i]);
		}
		//double[] specINT = spec.getSpecINT();
		// integer peptide mz
		long[] pepMZ = new long[pepMZDouble.length + 1];
	
		for (int i = 0; i < pepMZDouble.length; i++){
			pepMZ[i] = Math.round(pepMZDouble[i]);
		}
		pepMZ[pepMZDouble.length] = 0;
		// double peptide mz
		//ouble[] pepMZDouble = null;
		// intensity peptide double
		//double[] pepINT = null;
	   
		int matchedConsecutiveIon = 0;
		int matchedConsecutiveIonHighestValue = 0;
		int lastMatchedPeptideIndex = -1;
		
		//double[] newMatchedErrorList = null;
	
		double score = 0.0;
		// itType replacement --> an integer, spectrum MZ pointer
		// also replaces pType???
		int itType = 0;
		// m_lId replacement --> ist irrelevant
		// int m_lId = specMZ.length;
		// tTypeSize ist wieder ein max wert --> max length of spectrum
		int tTypeSize = specMZ.length;
		// tType replacement --> counter for spectrum
		int tType = 0;
		// tStep replacement --> spectrum peaks / peptid peaks
		int tStep = (int) (0.5 + ((double) ((double) specMZ.length / ((double) pepMZ.length))));
		if (tStep < 1) {
			tStep = 1;
		}
	
		// a ist ein wichtiger counter
		int a = 0;
		// b ist ein wichtiger counter
		int b = 0;
	
		// lcount ist die count im gepufferten scoring vektor
		int lCount = 0;
		// gepufferte intensit�ten nur f�r matches aus specINT und pepINT
		double[] m_pafI = new double[256];
		double[] m_pafSeq = new double[256];
	
		// kann boolean
		boolean lType = false;
	
		// HAUPTSCHLEIFE MATCHING
		while (pepMZ[a] != 0 && tType != tTypeSize) {
	
			// kann boolean
			lType = false;
	
			// ERSTES IF --> JUMPING
			if (specMZ[tType] < pepMZ[a]) {
				lType = true;
				while ((tType + tStep < tTypeSize) && (specMZ[tType + tStep] < pepMZ[a])) {
					tType += tStep;
					// pType += tStep;
				}
				do {
					tType++;
					// pType++;
					// KS: Solange weiter in einzelschritt bis �ber lSeq(PEPTID) pointer
				} while (tType < tTypeSize && specMZ[tType] < pepMZ[a]);
			} else if (specMZ[tType] > pepMZ[a]) {
				do {
					a++;
				} while (specMZ[tType] > pepMZ[a] && pepMZ[a] != 0);
			}
			// ZWEITES IF --> ABBRUCH
			if (pepMZ[a] == 0 || tType == tTypeSize) {
				break;
			}
			// DRITTES IF --> MATCH
			if (specMZ[tType] == pepMZ[a]) {
				m_pafI[lCount] = specINT[itType + tType];
				m_pafSeq[lCount] = pepINT[a];
	
				// code dealing with collection of errors
				//newMatchedErrorList[lCount] = pepMZDouble[a] - spec.getSpecMZ()[(int) (tType / 3)];
	
				// iterate lCount -> number of matches
				lCount++;
				// code dealing with consecutive matches
				if (a == lastMatchedPeptideIndex + 1) {
					matchedConsecutiveIon++;
				} else {
					if (matchedConsecutiveIon > matchedConsecutiveIonHighestValue) {
						matchedConsecutiveIonHighestValue = matchedConsecutiveIon;
					}
					matchedConsecutiveIon = 0;
				}
				lastMatchedPeptideIndex = a;
			}
			// VIERTES IF --> NACHBEREITUNG
			if (lType) {
				a++;
			} else {
				tType++;
			}
	
		}
	
		// SCORING
		double matchIntensity = 0.0;
		for (int i = 0; i < lCount; i++) {
			score += m_pafI[i] * m_pafSeq[i];
			matchIntensity += m_pafI[i];
		}
		
		Score result = new Score(score, lCount);
		// set values in PSM
		// Score is now dealt with in the score() function
		result.sumOfMatchedIntensities += matchIntensity;
		if (matchedConsecutiveIon > matchedConsecutiveIonHighestValue) {
			matchedConsecutiveIonHighestValue = matchedConsecutiveIon;
		}
	
		result.matchedIonIntensity = matchIntensity;
		result.matchedIonToTotal = matchIntensity / m_pafSeq.length;
		result.matchedIonNumLongestConsecutiveMatch = matchedConsecutiveIonHighestValue;
		//result.matchErrorListIon = newMatchedErrorList;

		return result;
	}

	private static double factorial(int number) {
		double f = 1.0;
		for (int i = 1; i <= number; i++) {
			f *= i;
		}
		return f;
	}

	public double calculate() {
		Future<Map<String, double[]>> ionsFuture = pool.submit(new PrettyFragmentationCallable(this.sequence));
		Map<String, double[]> result = null;
		try {
			result = ionsFuture.get();
		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		double[] y = result.get("y");
		double[] ypluplus = result.get("y++");
		double[] b = result.get("b");
		double[] bplusplus = result.get("b++");
		double[] theoreticalIntensity = new double[y.length];
		double[] theoreticalIntensity2 = new double[pepMZ.length];
		for (int i = 0; i < y.length; i++) {
			theoreticalIntensity[i] = 1.0;
		}
		for (int i = 0 ; i < this.pepMZ.length; i++){
			theoreticalIntensity2[i] = 1.0;
		}
		// scoring function call, separate for different ion types
		Score matched_B_1_ions = score(b, theoreticalIntensity, this.pepMZ, this.pepINT);
		Score matched_B_2_ions = score(bplusplus, theoreticalIntensity, this.pepMZ, this.pepINT);
		Score matched_Y_1_ions = score(y, theoreticalIntensity, this.pepMZ, this.pepINT);
		Score matched_Y_2_ions = score(ypluplus, theoreticalIntensity, this.pepMZ, this.pepINT);

		if (matched_B_1_ions.score != 0.0)
			return (Math.log10(factorial(matched_B_2_ions.matchedIons) * factorial(matched_Y_2_ions.matchedIons)
					* factorial(matched_B_1_ions.matchedIons) * factorial(matched_Y_1_ions.matchedIons) * matched_B_1_ions.score));
		return 0.0;
	}

    public double getAAMass(char c) {
        switch (c) {
            case 'A':
                return 71.03712;
            case 'B':
                return 114.1038; // This value is not used: see Note 1 at the end of this method.
            case 'C':
                return 103.00919;
            case 'c':
                return 103.00919 + 57.02146;
            case 'D':
                return 115.02695;
            case 'E':
                return 129.0426;
            case 'F':
                return 147.06842;
            case 'G':
                return 57.02147;
            case 'H':
                return 137.05891;
            case 'I':
                return 113.08407;
            case 'J':
                return 113.08407;
            case 'K':
                return 128.09497;
            case 'L':
                return 113.08407;
            case 'M':
                return 131.04049;
            case 'm':
                return 131.04049 + 15.99491;
            case 'N':
                return 114.04293;
            case 'O':
                return 237.31; // pyrolysine
            case 'P':
                return 97.05277;
            case 'Q':
                return 128.05858;
            case 'R':
                return 156.10112;
            case 'S':
                return 87.03203;
            case 'T':
                return 101.04768;
            case 'U':
                return 103.1388 - 32.066 + 78.96; // selenocysteine
            case 'V':
                return 99.06842;
            case 'W':
                return 186.07932;
            case 'X':
                return 0.0; // X is interpreted as an internal cleavage site/stop codon. It is the
                            // equivalent of "*".
            case 'Y':
                return 163.06333;
            case 'Z':
                return 128.1307;
        }
        return 0.0;
    }

}
