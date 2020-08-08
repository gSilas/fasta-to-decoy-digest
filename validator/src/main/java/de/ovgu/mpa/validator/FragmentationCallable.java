package de.ovgu.mpa.validator;

import java.util.Arrays;
import java.util.concurrent.Callable;

public class FragmentationCallable implements Callable<double[]> {

    private final double PROTON = 1.007276;

    private String pepSequence;

    public FragmentationCallable(String pepSequence) {
        this.pepSequence = pepSequence;
    }

    @Override
    public double[] call() throws Exception {
        
        double[] fragmentIons = new double[this.pepSequence.length() * 4];

        int a = 0;
        int lCount = 0;

        // DEAL WITH PROTEIN N-TERMINUS
        // dValue = 0.0;
        // deal with non-hydrolytic cleavage
        double bValue = 0.0;
        double yValue = 18.010560035000001;

        // MAIN LOOP
        while (a <= this.pepSequence.length() - 1) {
            bValue += FragmentationCalculator.getAAMass(this.pepSequence.charAt(a));
            yValue += FragmentationCalculator.getAAMass(this.pepSequence.charAt(this.pepSequence.length() - 1 - a));
            fragmentIons[lCount] = PROTON + bValue;
            lCount++;
            fragmentIons[lCount] = (PROTON * 2.0 + bValue) / 2.0;
            lCount++;
            fragmentIons[lCount] = PROTON + yValue;
            lCount++;
            fragmentIons[lCount] = (PROTON * 2.0 + yValue) / 2.0;
            lCount++;

            a++;
        }

        Arrays.sort(fragmentIons);
        return fragmentIons;
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
