package ru.spbau.mit.kazakov.Bioinf;

import javafx.util.Pair;
import org.apache.commons.io.FileUtils;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Bioinf {
    private static final String HCD_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.txt";
    private static final String CID_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.CID.txt";
    private static final String SCANS_DESC_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign";
    private static final String CHAINS_FILE = "Mab.fasta";
    private static final String HCD_ANALYSIS_FILE = "hcd";
    private static final String CID_ANALYSIS_FILE = "cid";
    private static final String RESOURCES_FOLDER = "src/main/resources/";

    private static int NUMBER_OF_HCD_TRIES = 2;
    private static int NUMBER_OF_CID_TRIES = 30;

    private static final Pattern activationPattern = Pattern.compile("ACTIVATION=(...)");
    private static final Pattern scanNumberPattern = Pattern.compile("SCANS=(\\d+)");
    private static final Pattern precursorMassPattern = Pattern.compile("PRECURSOR_MASS=(\\d+\\.\\d+)");
    private static final double epsilon = 0.0001;

    private static Map<Integer, Scan> scans = new HashMap<>();
    private static List<Pair<Scan, Scan>> pairs = new ArrayList<>();

    private static String heavyChain;
    private static String lightChain;


    public static void main(String[] args) throws IOException {
        readChainFile();

        readScanResFile(RESOURCES_FOLDER + HCD_SCANS_RES_FILE);
        readScanDescFile();
        printAnalysis(NUMBER_OF_CID_TRIES, NUMBER_OF_HCD_TRIES, HCD_ANALYSIS_FILE);

        scans.clear();
        pairs.clear();
        readScanResFile(RESOURCES_FOLDER + CID_SCANS_RES_FILE);
        readScanDescFile();
        printAnalysis(2, 2, CID_ANALYSIS_FILE);

    }

    private static void printAnalysis(int numberOfCidTries, int numberOfHcdTries, String file) throws FileNotFoundException {
        List<String> HCDrecognized = new ArrayList<>();
        List<String> CIDrecognized = new ArrayList<>();
        for (Scan scan : scans.values()) {
            int numberOfTries = scan.getActivation() == Activation.CID ? numberOfCidTries : numberOfHcdTries;
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

        try (PrintWriter hcdWriter = new PrintWriter(RESOURCES_FOLDER + file)) {
            hcdWriter.println("HCD:");
            printList(hcdWriter, HCDrecognized);

            hcdWriter.println();
            hcdWriter.println("CID:");
            printList(hcdWriter, CIDrecognized);

            hcdWriter.println();
            printPairs(hcdWriter);
        }
    }

    private static void printPairs(PrintWriter out) {
        out.println("Pairs:");
        for(Pair<Scan, Scan> pair : pairs) {
            if (!pair.getKey().getSequences().isEmpty() && !pair.getValue().getSequences().isEmpty()) {
                String sequence = pair.getKey().getSequences().get(0).getSequence();
                if (sequence.equals(pair.getValue().getSequences().get(0).getSequence())) {
                    out.println();
                    out.println("SCAN_1: " + pair.getKey().getNumber());
                    out.println("SCAN_2: " + pair.getValue().getNumber());
                    out.println("Sequence: " + sequence);
                    out.println("RnkScr_1: " + pair.getKey().getSequences().get(0).getRnkScr());
                    out.println("RnkScr_1: " + pair.getValue().getSequences().get(0).getRnkScr());
                }
            }
        }
    }
    
    private static void printList(PrintWriter out, List<String> list) {
        for (String string : list) {
            out.println(string);
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

    @NotNull
    private static Sequence parseSequenceDesc(@NotNull String sequenceDesc) {
        String[] desc = sequenceDesc.split("\\s");
        return new Sequence(desc[7], Double.parseDouble(desc[1]));
    }

    private static void parseScanDesc(@NotNull String scanDesc) {
        String[] parameters = scanDesc.trim().split("\\n");

        Matcher numMatcher = scanNumberPattern.matcher(parameters[2]);
        numMatcher.find();
        int num = Integer.parseInt(numMatcher.group(1));

        Matcher activationMatcher = activationPattern.matcher(parameters[3]);
        activationMatcher.find();
        String activation = activationMatcher.group(1);

        Matcher massMatcher = precursorMassPattern.matcher(parameters[6]);
        massMatcher.find();
        double mass = Double.parseDouble(massMatcher.group(1));

        scans.get(num).setActivation(Activation.getActivationByString(activation));
        scans.get(num).setPrecursorMass(mass);

        Scan pair = scans.get(num - 1);
        if (pair != null && Math.abs(pair.getPrecursorMass() - mass) <= epsilon) {
            pairs.add(new Pair<>(pair, scans.get(num)));
        }
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

