package de.ovgu.mpa.validator.fasta;

import de.ovgu.mpa.validator.ValidatorConfig;

public class CompareResult {

    // <length, target proteins>
    public long[] lengthTargetArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];
    // <length, decoy proteins>
    public long[] lengthDecoyArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];
    // <length, bins>
    public long[][] lengthBinMatrix = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH][11];
    // <length, great matches>
    // HashMap<Integer, Integer> lengthMap = new HashMap<>();
    public long[] lengthMatchArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];
    
    public long[] cosineSimilarityBins = new long[11];

    public long matchCounter;
    public long greatMatchCounter;

    public long[] bestCosineSimilarityBins = new long[11];
    public long[][] lengthBestBinMatrix = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH][11];

    public long[] fragmentMatch = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH * 4];

    public long[] lengthDuplicateArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];

}
