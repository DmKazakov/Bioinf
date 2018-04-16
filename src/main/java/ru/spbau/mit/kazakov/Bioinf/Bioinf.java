package ru.spbau.mit.kazakov.Bioinf;

import org.apache.commons.io.FileUtils;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class Bioinf {
    private static final String HCD_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.txt";
    private static final String CID_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.CID.txt";
    private static final String SCANS_DESC_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign";
    private static final String CHAINS_FILE = "Mab.fasta";
    private static final String RESOURCES_FOLDER = "src/main/resources/";
    private static int NUMBER_OF_HCD_TRIES = 2;
    private static int NUMBER_OF_CID_TRIES = 30;
    private static Map<Integer, Scan> scans = new HashMap<>();

    private static String heavyChain;
    private static String lightChain;


    public static void main(String[] args) throws IOException {
        readScanResFile(RESOURCES_FOLDER + HCD_SCANS_RES_FILE);
        readScanDescFile();
        readChainFile();

        List<String> HCDrecognized = new ArrayList<>();
        List<String> CIDrecognized = new ArrayList<>();
        for (Scan scan : scans.values()) {
            int numberOfTries = scan.getActivation() == Activation.CID ? NUMBER_OF_CID_TRIES : NUMBER_OF_HCD_TRIES;
            for (int i = 0; i < Math.min(numberOfTries + 1, scan.getSequences().size()); i++) {
                Sequence sequence = scan.getSequences().get(i);
                if (isSubstring(lightChain, sequence.getSequence()) || isSubstring(heavyChain, sequence.getSequence())) {
                    if (scan.getActivation().equals(Activation.HCD)) {
                        HCDrecognized.add(sequence.getSequence());
                    } else {
                        CIDrecognized.add(sequence.getSequence());
                    }
                }
            }
        }

        System.out.println("HCD:");
        for (String sequence : HCDrecognized) {
            System.out.println(sequence);
        }
        System.out.println("\nCID:");
        for (String sequence : CIDrecognized) {
            System.out.println(sequence);
        }
    }

    private static void readScanDescFile() throws IOException {
        File scansDescFile = new File(RESOURCES_FOLDER + SCANS_DESC_FILE);
        String scansDescFileContent = FileUtils.readFileToString(scansDescFile).trim();
        String[] scansDesc = scansDescFileContent.split("END IONS");
        for (String scanDesc : scansDesc) {
            parseScanDesc(scanDesc.trim());
        }
    }

    private static void readChainFile() throws IOException {
        File chainsFile = new File(RESOURCES_FOLDER + CHAINS_FILE);
        String[] chains = FileUtils.readFileToString(chainsFile).trim().split("\\n\\n");
        List<String> lightChainLines = Arrays.asList(chains[0].split("\n"));
        lightChain = String.join("", lightChainLines.subList(1, lightChainLines.size()));
        List<String> heavyChainLines = Arrays.asList(chains[1].split("\n"));
        heavyChain = String.join("", heavyChainLines.subList(1, heavyChainLines.size()));
    }

    private static void readScanResFile(String path) throws IOException {
        File scanResFile = new File(path);
        String scanResFileContent = FileUtils.readFileToString(scanResFile).trim();
        String[] scansRes = scanResFileContent.split("\\n\\n");
        for (String scanRes : scansRes) {
            parseScanRes(scanRes.trim());
        }
    }

    private static void parseScanRes(@NotNull String scanRes) {
        String[] sequencesDesc = scanRes.split("\\n");
        if (sequencesDesc[0].charAt(0) != '>') {
            return;
        }
        if (!sequencesDesc[0].split("\\s")[2].matches("\\d+")) {
            for (String desc : sequencesDesc) {
                if (desc.charAt(0) == '#') {
                    continue;
                }
                String stringNum = desc.split("\\s")[2];
                int num = Integer.parseInt(stringNum.substring(0, stringNum.length() - 8));
                scans.put(num, new Scan(num, new ArrayList<>()));
            }
            return;
        }

        int num = Integer.parseInt(sequencesDesc[0].split("\\s")[2]);
        if (sequencesDesc[1].equals("# No solutions found.")) {
            scans.put(num, new Scan(num, new ArrayList<>()));
        } else {
            List<Sequence> sequences = new ArrayList<>();
            for (String sequenceDesc : Arrays.asList(sequencesDesc).subList(2, sequencesDesc.length)) {
                sequences.add(parseSequenceDesc(sequenceDesc));
            }
            scans.put(num, new Scan(num, sequences));
        }
    }

    private static Sequence parseSequenceDesc(@NotNull String sequenceDesc) {
        String[] desc = sequenceDesc.split("\\s");
        return new Sequence(desc[7], Double.parseDouble(desc[1]));
    }

    private static void parseScanDesc(@NotNull String scanDesc) {
        String[] parameters = scanDesc.trim().split("\\n");
        int num = Integer.parseInt(parameters[2].substring(6, parameters[2].length() - 1));
        String activation = parameters[3].substring(11, parameters[3].length() - 1);

        scans.get(num).setActivation(Activation.getActivationByString(activation));
    }

    private static boolean isSubstring(@NotNull String string, @NotNull String substring) {
        String concat = substring + "#" + string;
        int[] prefixFunction = new int[concat.length()];

        for (int i = 1; i < concat.length(); i++) {
            int count = prefixFunction[i - 1];
            while (count > 0 && concat.charAt(count) != concat.charAt(i))
                count = prefixFunction[count - 1];
            if (concat.charAt(count) == concat.charAt(i))
                count++;
            prefixFunction[i] = count;
        }

        for (int i = 0; i < concat.length(); i++) {
            if (prefixFunction[i] == substring.length()) {
                return true;
            }
        }

        return false;
    }
}

