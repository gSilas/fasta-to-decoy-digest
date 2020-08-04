package de.ovgu.mpa.validator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class CSVWriter {
    
    public CSVWriter() {}

    public String convertToCSV(String[] data) {
        return Stream.of(data)
          .map(this::escapeSpecialCharacters)
          .collect(Collectors.joining(","));
    }

    public void createCSV(List<String[]> dataLines, String csvName) throws FileNotFoundException {
        File csvOutputFile = new File(csvName);
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            dataLines.stream()
            .map(this::convertToCSV)
            .forEach(pw::println);
        }
        assert(csvOutputFile.exists());
    }

    public String escapeSpecialCharacters(String data) {
        String escapedData = data.replaceAll("\\R", " ");
        if (data.contains(",") || data.contains("\"") || data.contains("'")) {
            data = data.replace("\"", "\"\"");
            escapedData = "\"" + data + "\"";
        }
        return escapedData;
    }
}