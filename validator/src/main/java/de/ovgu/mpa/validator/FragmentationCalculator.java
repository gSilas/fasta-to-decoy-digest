package de.ovgu.mpa.validator;

import java.util.Arrays;

public class FragmentationCalculator {

	private static double m_dProton = 1.007276;
	
    public enum IonSeries {
		B_ION_1,
		B_ION_2,
		Y_ION_1,
		Y_ION_2;
	}

    public enum AminoAcidMasses {
		A('A', 71.03712),
		B('B', 114.1038),	// This value is not used: see Note 1 at the end of this method.
		C('C', 103.00919),
		
		C_CAM('c', 103.00919 + 57.02146),
		
		D('D', 115.02695),
		E('E', 129.0426),
		F('F', 147.06842),
		G('G', 57.02147),
		H('H', 137.05891),
		I('I', 113.08407),
		J('J', 113.08407),
		K('K', 128.09497),
		L('L', 113.08407),
		M('M', 131.04049),
		
		M_O('m', 131.04049 + 15.99491),
		
		N('N', 114.04293),
		O('O', 237.31),	// pyrolysine
		P('P', 97.05277),
		Q('Q', 128.05858),
		R('R', 156.10112),
		S('S', 87.03203),
		T('T', 101.04768),
		U('U', 103.1388 - 32.066 + 78.96),	// selenocysteine
		V('V', 99.06842),
		W('W', 186.07932),
		X('X', 0.0),	// X is interpreted as an internal cleavage site/stop codon. It is the equivalent of "*".
		Y('Y', 163.06333),
		Z('Z', 128.1307); // This value is not used: see Note 1 at the end of this method
		
		private char oneLetter;
		private double mass;

		private AminoAcidMasses(char o, double mass) {
			this.oneLetter = o;
			this.mass = mass;
		}

		public static double massFromOneLetter(char c) {
			for (AminoAcidMasses amac : AminoAcidMasses.values()) {
				if (amac.oneLetter == c) {
					return amac.mass;
				}
			}
			return 0.0;
		}

	}
	
	public static double getAAMass(char c) {
		switch (c) {
			case 'A': return 71.03712;
			case 'B': return 114.1038;	// This value is not used: see Note 1 at the end of this method.
			case 'C': return 103.00919;
			case 'c': return 103.00919 + 57.02146;
			case 'D': return 115.02695;
			case 'E': return 129.0426;
			case 'F': return 147.06842;
			case 'G': return 57.02147;
			case 'H': return 137.05891;
			case 'I': return 113.08407;
			case 'J': return 113.08407;
			case 'K': return 128.09497;
			case 'L': return 113.08407;
			case 'M': return 131.04049;
			case 'm': return 131.04049 + 15.99491;
			case 'N': return 114.04293;
			case 'O': return 237.31;	// pyrolysine
			case 'P': return 97.05277;
			case 'Q': return 128.05858;
			case 'R': return 156.10112;
			case 'S': return 87.03203;
			case 'T': return 101.04768;
			case 'U': return 103.1388 - 32.066 + 78.96;	// selenocysteine
			case 'V': return 99.06842;
			case 'W': return 186.07932;
			case 'X': return 0.0;	// X is interpreted as an internal cleavage site/stop codon. It is the equivalent of "*".
			case 'Y': return 163.06333;
			case 'Z': return 128.1307; 
		}
		return 0.0;
	}
    
    public static void getFragmentIons(String peptideSequence, double[] fragmentIons) {
		//int startIndex = 0;
		//startIndex = add_B(peptideSequence, fragmentIons, 1.0, startIndex);
		//startIndex = add_B(peptideSequence, fragmentIons, 2.0, startIndex);
		//startIndex = add_Y(peptideSequence, fragmentIons, 1.0, startIndex);
		//startIndex = add_Y(peptideSequence, fragmentIons, 2.0, startIndex);
		fragmentIonBY(peptideSequence, fragmentIons);
		Arrays.sort(fragmentIons);
	}
	

	private static void fragmentIonBY(String pepSeq, double[] fragmentIons) {
		int a = 0;
		int lCount = 0;
		
		// DEAL WITH PROTEIN N-TERMINUS
		// dValue = 0.0;
		// deal with non-hydrolytic cleavage
		double bValue = 0.0;
		double yValue = 18.010560035000001;
		
		// MAIN LOOP
		while (a <= pepSeq.length() - 1) {
			bValue += FragmentationCalculator.getAAMass(pepSeq.charAt(a));
			yValue += FragmentationCalculator.getAAMass(pepSeq.charAt(pepSeq.length() -1 - a));
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(bValue, 1.0);
			lCount++;
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(bValue, 2.0);
			lCount++;
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(yValue, 1.0);
			lCount++;
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(yValue, 2.0);
			lCount++;

			a++;
		}

		while (lCount < ValidatorConfig.MAXIMUM_PEP_LENGTH * 4) {
			fragmentIons[lCount] = Double.POSITIVE_INFINITY;
			lCount++;
		}

	}
    
    private static int add_B(String pepSeq, double[] fragmentIons, double charge, int startIndex) {
		// inits
		int a = 0;
		// double dValue = 036;
		int lCount = startIndex;
		
		// DEAL WITH PROTEIN N-TERMINUS
		// dValue = 0.0;
		// deal with non-hydrolytic cleavage
		double dValue = 0.0;
		
		// MAIN LOOP
		while (a <= pepSeq.length() - 1) {
			// add amino acid mass
			dValue += AminoAcidMasses.massFromOneLetter(pepSeq.charAt(a));
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(dValue, charge);
//			dValue -= m_dProton;
			// TODO: DEAL WITH MODIFICATIONS
			// convert to integer
			lCount++;
			a++;
		}
		
		return lCount;
	}

	private static int add_Y(String pepSeq, double[] fragmentIons, double charge, int startIndex) {
		// initial stuff C6H12ON2
		// inits
		int a = pepSeq.length() - 1;
		// double dValue = 18.0105647;
		boolean bZero = false;
		int lCount = startIndex;
		
		// deal with non-hydrolytic cleavage
		// dValue = 18.010560035000001;
		// deal with protein C-teminus
		double dValue = 18.010560035000001;

		// MAIN LOOP
		while (a >= 0) {
			// add amino acid mass
			dValue += AminoAcidMasses.massFromOneLetter(pepSeq.charAt(a));
			fragmentIons[lCount] = FragmentationCalculator.mconvert_double(dValue, charge);
			if (bZero)	{
				if(a < 5)	{
					lCount++;
				}
			}
			else	{
				lCount++;
			}
			// end of loop down iteration
			a--;
		}
		return lCount;
    }
    
    private static double mconvert_double(double mass, double charge) {
		return (m_dProton*charge + mass) / charge;
	}
	
}