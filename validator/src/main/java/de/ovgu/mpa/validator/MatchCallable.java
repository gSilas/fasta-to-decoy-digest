package de.ovgu.mpa.validator;

import java.util.concurrent.Callable;
import java.util.concurrent.Future;

public class MatchCallable implements Callable<Double> {

    Future<double[]> decoyIonsFuture;
    Future<double[]> targetIonsFuture;
    double tolerance;

    public MatchCallable(Future<double[]> decoyIonsFuture, Future<double[]> targetIonsFuture, double tolerance) {
        this.decoyIonsFuture = decoyIonsFuture;
        this.targetIonsFuture = targetIonsFuture;
        this.tolerance = tolerance;
    }

    @Override
    public Double call() throws Exception {

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

        double cosineSimilarity = fragmentMatch / (Math.sqrt(targetIons.length) * Math.sqrt(decoyIons.length));

        return cosineSimilarity;
    }
}
