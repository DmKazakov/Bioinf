package ru.spbau.mit.kazakov.Bioinf;

import javafx.util.Pair;
import org.apache.commons.collections4.Bag;
import org.apache.commons.collections4.bag.TreeBag;
import org.apache.commons.io.FileUtils;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class Bioinf {
    private static final String SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.txt";
    private static final String CID_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.CID.txt";
    private static final String HCD_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.HCD.txt";
    private static final String SCANS_DESC_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign";
    private static final String CHAINS_FILE = "Mab.fasta";
    private static final String FIRST_ANALYSIS_FILE = "1";
    private static final String SECOND_ANALYSIS_FILE = "2";
    private static final String THIRD_ANALYSIS_FILE = "3";
    private static final String RESOURCES_FOLDER = "src/main/resources/";

    private static int ALL = 20;

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

        readScanResFile(RESOURCES_FOLDER + SCANS_RES_FILE);
        readScanDescFile();
        Pair<Coverage, Coverage> firstFile = printStat(FIRST_ANALYSIS_FILE);

        scans.clear();
        pairs.clear();
        readScanResFile(RESOURCES_FOLDER + CID_SCANS_RES_FILE);
        readScanDescFile();
        Pair<Coverage, Coverage> secondFile = printStat(SECOND_ANALYSIS_FILE);

        scans.clear();
        pairs.clear();
        readScanResFile(RESOURCES_FOLDER + HCD_SCANS_RES_FILE);
        readScanDescFile();
        Pair<Coverage, Coverage> thirdFile = printStat(THIRD_ANALYSIS_FILE);

        System.out.println("1 file, HCD 20: " + firstFile.getValue().getHcdHeavyChain());
        System.out.println("1 file, HCD 02: " + firstFile.getKey().getHcdHeavyChain());
        System.out.println("2 file, HCD 20: " + secondFile.getValue().getHcdHeavyChain());
        System.out.println("2 file, HCD 02: " + secondFile.getKey().getHcdHeavyChain());
        System.out.println("3 file, HCD 20: " + thirdFile.getValue().getHcdHeavyChain());
        System.out.println("3 file, HCD 02: " + thirdFile.getKey().getHcdHeavyChain());
        System.out.println("1 file, CID 20: " + firstFile.getValue().getCidHeavyChain());
        System.out.println("1 file, CID 02: " + firstFile.getKey().getCidHeavyChain());
        System.out.println("2 file, CID 20: " + secondFile.getValue().getCidHeavyChain());
        System.out.println("2 file, CID 02: " + secondFile.getKey().getCidHeavyChain());
        System.out.println("3 file, CID 20: " + thirdFile.getValue().getCidHeavyChain());
        System.out.println("3 file, CID 02: " + thirdFile.getKey().getCidHeavyChain());
        System.out.println();

        System.out.println("1 file, HCD 20: " + firstFile.getValue().getHcdLightChain());
        System.out.println("1 file, HCD 02: " + firstFile.getKey().getHcdLightChain());
        System.out.println("2 file, HCD 20: " + secondFile.getValue().getHcdLightChain());
        System.out.println("2 file, HCD 02: " + secondFile.getKey().getHcdLightChain());
        System.out.println("3 file, HCD 20: " + thirdFile.getValue().getHcdLightChain());
        System.out.println("3 file, HCD 02: " + thirdFile.getKey().getHcdLightChain());
        System.out.println("1 file, CID 20: " + firstFile.getValue().getCidLightChain());
        System.out.println("1 file, CID 02: " + firstFile.getKey().getCidLightChain());
        System.out.println("2 file, CID 20: " + secondFile.getValue().getCidLightChain());
        System.out.println("2 file, CID 02: " + secondFile.getKey().getCidLightChain());
        System.out.println("3 file, CID 20: " + thirdFile.getValue().getCidLightChain());
        System.out.println("3 file, CID 02: " + thirdFile.getKey().getCidLightChain());
    }

    @NotNull
    private static Pair<Coverage, Coverage> printStat(@NotNull String file) throws FileNotFoundException {
        Coverage coverageTwoInterpretations;
        Coverage coverageAllInterpretations;

        try (PrintWriter writer = new PrintWriter(RESOURCES_FOLDER + file)) {
            coverageTwoInterpretations = printMatches(writer, 2, 2);
            writer.println();
            coverageAllInterpretations = printMatches(writer, ALL, ALL);
            writer.println();
            printFirstInterpretationStat(writer, 1);
            writer.println();
            writer.println();
            printPairs(writer, 1);
            writer.println();
            printPairs(writer, 2);
            writer.println();
            printPairs(writer, ALL);
        }

        return new Pair<>(coverageTwoInterpretations, coverageAllInterpretations);
    }

    private static void printFirstInterpretationStat(PrintWriter out, int numberOfInterpretations) {
        int hcdCorrectInterpretations = 0;
        int cidCorrectInterpretations = 0;
        int hcdScans = 0;
        int cidScans = 0;
        for (Scan scan : scans.values()) {
            if (scan.getActivation().equals(Activation.CID)) {
                cidScans++;
            } else {
                hcdScans++;
            }
            for (int i = 0; i < Math.min(numberOfInterpretations + 1, scan.getSequences().size()); i++) {
                Sequence sequence = scan.getSequences().get(i);
                if (!isSubstring(lightChain, sequence.getSequence()).isEmpty()) {
                    if (scan.getActivation().equals(Activation.HCD)) {
                        hcdCorrectInterpretations++;
                    } else {
                        cidCorrectInterpretations++;
                    }
                } else if (!isSubstring(heavyChain, sequence.getSequence()).isEmpty()) {
                    if (scan.getActivation().equals(Activation.HCD)) {
                        hcdCorrectInterpretations++;
                    } else {
                        cidCorrectInterpretations++;
                    }
                }

            }
        }
        out.println("HCD: " + ((double) Math.round(1000d * (hcdCorrectInterpretations / (double) hcdScans)) / 10d)
                + "% of correct first interpretations");
        out.println("CID: " + ((double) Math.round(1000d * (cidCorrectInterpretations / (double) cidScans)) / 10d)
                + "% of correct first interpretations");
    }

    private static Coverage printMatches(PrintWriter out, int numberOfCidTries, int numberOfHcdTries) {
        Bag<Pair<Integer, String>> hcdLightChain = new TreeBag<>(Comparator.comparing(Pair::getKey));
        Bag<Pair<Integer, String>> hcdHeavyChain = new TreeBag<>(Comparator.comparing(Pair::getKey));
        Bag<Pair<Integer, String>> cidLightChain = new TreeBag<>(Comparator.comparing(Pair::getKey));
        Bag<Pair<Integer, String>> cidHeavyChain = new TreeBag<>(Comparator.comparing(Pair::getKey));

        for (Scan scan : scans.values()) {
            int numberOfTries = scan.getActivation() == Activation.CID ? numberOfCidTries : numberOfHcdTries;
            for (int i = 0; i < Math.min(numberOfTries + 1, scan.getSequences().size()); i++) {
                Sequence sequence = scan.getSequences().get(i);
                List<Integer> positions = isSubstring(lightChain, sequence.getSequence());
                if (!positions.isEmpty()) {
                    if (scan.getActivation().equals(Activation.HCD)) {
                        for (int position : positions) {
                            hcdLightChain.add(new Pair<>(position, sequence.getSequence()));
                        }
                    } else {
                        for (int position : positions) {
                            cidLightChain.add(new Pair<>(position, sequence.getSequence()));
                        }
                    }
                }
                positions = isSubstring(heavyChain, sequence.getSequence());
                if (!positions.isEmpty()) {
                    if (scan.getActivation().equals(Activation.HCD)) {
                        for (int position : positions) {
                            hcdHeavyChain.add(new Pair<>(position, sequence.getSequence()));
                        }
                    } else {
                        for (int position : positions) {
                            cidHeavyChain.add(new Pair<>(position, sequence.getSequence()));
                        }
                    }
                }
            }
        }

        BitSet coverage = getCoverage(hcdLightChain, lightChain.length());
        String hcdLightChainCoverage = highlightCoverage(coverage, lightChain);
        out.println("HCD light chain (" + numberOfHcdTries + " interpretations, "
                + countCoverage(coverage, lightChain.length()) + "% coverage):");
        out.println(lightChain);
        printFoundMatches(out, hcdLightChain);
        out.println();

        coverage = getCoverage(hcdHeavyChain, heavyChain.length());
        String hcdHeavyChainCoverage = highlightCoverage(coverage, heavyChain);
        out.println("HCD heavy chain (" + numberOfHcdTries + " interpretations, "
                + countCoverage(coverage, heavyChain.length()) + "% coverage):");
        out.println(heavyChain);
        printFoundMatches(out, hcdHeavyChain);
        out.println();

        coverage = getCoverage(cidLightChain, lightChain.length());
        String cidLightChainCoverage = highlightCoverage(coverage, lightChain);
        out.println("CID light chain (" + numberOfHcdTries + " interpretations, "
                + countCoverage(coverage, lightChain.length()) + "% coverage):");
        out.println(lightChain);
        printFoundMatches(out, cidLightChain);
        out.println();

        coverage = getCoverage(cidLightChain, heavyChain.length());
        String cidHeavyChainCoverage = highlightCoverage(coverage, heavyChain);
        out.println("CID heavy chain (" + numberOfHcdTries + " interpretations, "
                + countCoverage(coverage, heavyChain.length()) + "% coverage):");
        out.println(heavyChain);
        printFoundMatches(out, cidHeavyChain);
        out.println();

        return new Coverage(hcdHeavyChainCoverage, cidHeavyChainCoverage, hcdLightChainCoverage, cidLightChainCoverage);
    }

    private static String highlightCoverage(BitSet coverage, String lightChain) {
        StringBuilder highlighted = new StringBuilder();

        for (int i = 0; i < lightChain.length(); i++) {
            if(coverage.get(i)) {
                highlighted.append(ConsoleColors.GREEN).append(lightChain.charAt(i)).append(ConsoleColors.RESET);
            } else {
                highlighted.append(lightChain.charAt(i));
            }
        }

        return highlighted.toString();
    }

    private static void printFoundMatches(@NotNull PrintWriter out, @NotNull Bag<Pair<Integer, String>> matches) {
        for (Pair<Integer, String> match : matches.uniqueSet()) {
            printWhitespaces(out, match.getKey());
            out.print(match.getValue());

            int count = matches.getCount(match);
            if (count > 1) {
                out.print("(" + count + ")");
            }

            out.println();
        }
    }

    private static void printWhitespaces(@NotNull PrintWriter out, int num) {
        for (int i = 0; i < num; i++) {
            out.print(' ');
        }
    }

    private static void printPairs(@NotNull PrintWriter out, int numberOfTries) {
        out.println(numberOfTries + " Pairs:");
        for (Pair<Scan, Scan> pair : pairs) {
            if (!pair.getKey().getSequences().isEmpty() && !pair.getValue().getSequences().isEmpty()) {
                IndexedSet<String> sequences = pair.getKey().getSequences().stream().limit(numberOfTries)
                        .map(Sequence::getSequence).collect(Collectors.toCollection(IndexedSet::new));
                List<String> itSequences = pair.getValue().getSequences().stream().limit(numberOfTries).map(Sequence::getSequence).collect(Collectors.toList());
                for (int i = 0; i < itSequences.size(); i++) {
                    if (sequences.getIndex(itSequences.get(i)) >= 0) {
                        out.println("SCAN_1: " + pair.getKey().getNumber());
                        out.println("SCAN_2: " + pair.getValue().getNumber());
                        out.println("Sequence: " + itSequences.get(i));
                        out.println("RnkScr_1: " + pair.getKey().getSequences().get(0).getRnkScr());
                        out.println("RnkScr_2: " + pair.getValue().getSequences().get(0).getRnkScr());
                        out.println("Index_1: " + sequences.getIndex(itSequences.get(i)));
                        out.println("Index_2: " + i);
                        out.println();
                    }
                }
            }
        }
    }

    private static void printList(@NotNull PrintWriter out, @NotNull List<String> list) {
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

    private static void readScanResFile(@NotNull String path) throws IOException {
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

    @NotNull
    private static List<Integer> isSubstring(@NotNull String string, @NotNull String substring) {
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

        List<Integer> ans = new ArrayList<>();
        for (int i = 0; i < concat.length(); i++) {
            if (prefixFunction[i] == substring.length()) {
                ans.add(i - 2 * substring.length());
            }
        }

        return ans;
    }

    private static BitSet getCoverage(@NotNull Bag<Pair<Integer, String>> matches, int stringLength) {
        BitSet covered = new BitSet(stringLength);

        for (Pair<Integer, String> match : matches) {
            covered.set(match.getKey(), match.getKey() + match.getValue().length());
        }

        return covered;
    }

    private static double countCoverage(@NotNull BitSet covered, int stringLength) {
        return (double) Math.round(1000d * (covered.cardinality() / (double) stringLength)) / 10d;
    }
}

