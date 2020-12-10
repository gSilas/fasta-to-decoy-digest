package de.ovgu.mpa.validator.mascot;

import java.util.HashMap;
import java.util.Map;

public class MascotCompareResult {

    // <length, score, match >
    // long[] lengthScoreArray = new long[ValidatorConfig.MAXIMUM_PEP_LENGTH];

    public Map<String, Double> mascotScoreMap = new HashMap<String, Double>();
    public Map<String, Double> cosineSimilarityMap = new HashMap<String, Double>();
    public Map<String, Integer> fragmentMatchMap = new HashMap<String, Integer>();
    public Map<String, Double> xTandemMap = new HashMap<String, Double>();

}
