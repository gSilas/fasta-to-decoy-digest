package de.ovgu.mpa.validator;

import static java.nio.file.StandardCopyOption.*;
import org.apache.commons.cli.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class App {

    public static void main(String[] args) {

        Options options = new Options();

        Option fasta = new Option("f", "read_fasta", true, "fasta file path");
        fasta.setRequired(false);
        options.addOption(fasta);

        Option compare = new Option("c", "compare_dbs", true, "paths to db1 abd db2");
        compare.setRequired(false);
        options.addOption(compare);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("Fasta DB Compare", options);
            System.exit(1);
        }

        if (cmd.hasOption("read_fasta")) {
            File batchDir = new File("tmp_batches");
            if (!batchDir.exists())
                batchDir.mkdir();
            File fastaFolder = new File("fasta");
            if (!fastaFolder.exists())
                fastaFolder.mkdir();

            // Create Target and Decoy Folder
            String fileName = cmd.getOptionValue("read_fasta")
                    .split("/")[cmd.getOptionValue("read_fasta").split("/").length - 1].split("\\.")[0];
            File targetFolder = new File(fastaFolder.getPath() + "/" + fileName);
            if (!targetFolder.exists())
                targetFolder.mkdir();
            File decoyFolder = new File(fastaFolder.getPath() + "/" + fileName + "_decoy");
            if (!decoyFolder.exists())
                decoyFolder.mkdir();

            processFasta(targetFolder.toString(), decoyFolder.toString(), cmd.getOptionValue("read_fasta"),
                    batchDir.toString(), fastaFolder);
        } else if (cmd.hasOption("compare_dbs")) {
            File resultsFolder = new File("results");
            if (!resultsFolder.exists())
                resultsFolder.mkdir();

            String targetFolder = cmd.getOptionValue("compare_dbs").split(" ")[0];
            String decoyFolder = cmd.getOptionValue("compare_dbs").split(" ")[1];
            compareDB(targetFolder.toString(), decoyFolder.toString(), resultsFolder.toString());

        } else {
            formatter.printHelp("Fasta DB Compare", options);
            System.exit(1);
        }
    }

    public static void processFasta(String targetFolder, String decoyFolder, String fastaPath, String batchFolder,
            File fastaFolder) {
        // Copy Fasta
        Path source = Paths.get(fastaPath);
        Path target = Paths.get(targetFolder + "/" + fastaPath.split("/")[fastaPath.split("/").length - 1]);
        try {
            Files.copy(source, target, REPLACE_EXISTING);
        } catch (IOException e2) {
            e2.printStackTrace();
        }
        // process fasta
        PeptideWriter writer = new PeptideWriter();
        try {
            writer.createFiles(targetFolder, batchFolder, target.toFile());
        } catch (IOException e) {
            e.printStackTrace();
        }

        // create decoy
        DecoyGenerator.generateReverseDecoyDatabase(target,
                Paths.get(decoyFolder + "/decoy_" + fastaPath.split("/")[fastaPath.split("/").length - 1]));

        writer = new PeptideWriter();
        File batchDir = new File("tmp_batches");
        if (!batchDir.exists())
            batchDir.mkdir();

        try {
            writer.createFiles(decoyFolder, batchFolder, Paths.get(decoyFolder + "/decoy_" + fastaPath.split("/")[fastaPath.split("/").length - 1]).toFile());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void compareDB(String targetFolder, String decoyFolder, String resultsFolder) {
        Statistics stats = new Statistics();
        try {
            stats.comparePeptides(decoyFolder, targetFolder, resultsFolder);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
