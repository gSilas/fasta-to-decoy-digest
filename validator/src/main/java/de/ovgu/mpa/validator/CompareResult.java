package de.ovgu.mpa.validator;

public class CompareResult {

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

    long matchCounter;
    long greatMatchCounter;

    long[] bestCosineSimilarityBins = new long[11];
    long[][] lengthBestBinMatrix = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH][11];

    long[] fragmentMatch = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];

}
