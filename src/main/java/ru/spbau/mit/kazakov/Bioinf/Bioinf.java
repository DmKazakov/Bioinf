package ru.spbau.mit.kazakov.Bioinf;

import com.google.common.collect.Sets;
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

import static ru.spbau.mit.kazakov.Bioinf.Activation.CID;
import static ru.spbau.mit.kazakov.Bioinf.Activation.HCD;

public class Bioinf {
    private static final String SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.txt";
    private static final String CID_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.CID.txt";
    private static final String HCD_SCANS_RES_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.HCD.txt";
    private static final String SCANS_DESC_FILE = "120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign";
    private static final String CHAINS_FILE = "Mab.fasta";
    private static final String FIRST_ANALYSIS_FILE = "1";
    private static final String SECOND_ANALYSIS_FILE = "2";
    private static final String THIRD_ANALYSIS_FILE = "3";
    private static final String GENERAL_STAT_FILE = "general";
    private static final String RESOURCES_FOLDER = "src/main/resources/";

    private static int ALL = 20;
    private static int DIFF = 10;

    private static final Pattern activationPattern = Pattern.compile("ACTIVATION=(...)");
    private static final Pattern scanNumberPattern = Pattern.compile("SCANS=(\\d+)");
    private static final Pattern precursorMassPattern = Pattern.compile("PRECURSOR_MASS=(\\d+\\.\\d+)");
    private static final double EPSILON = 0.0001;

    private static Map<Integer, Scan> scans = new HashMap<>();
    private static List<Pair<Scan, Scan>> pairs = new ArrayList<>();
    private static List<Set<Scan>> peptide = new ArrayList<>();

    private static String heavyChain;
    private static String lightChain;


    public static void main(String[] args) throws IOException {
        readChainFile();

        readScanResFile(RESOURCES_FOLDER + SCANS_RES_FILE);
        readScanDescFile();
        savePeptide();
        Pair<Coverage, Coverage> firstFile = printStat(FIRST_ANALYSIS_FILE);

        scans.clear();
        pairs.clear();
        readScanResFile(RESOURCES_FOLDER + CID_SCANS_RES_FILE);
        readScanDescFile();
        savePeptide();
        Pair<Coverage, Coverage> secondFile = printStat(SECOND_ANALYSIS_FILE);

        scans.clear();
        pairs.clear();
        readScanResFile(RESOURCES_FOLDER + HCD_SCANS_RES_FILE);
        readScanDescFile();
        savePeptide();
        Pair<Coverage, Coverage> thirdFile = printStat(THIRD_ANALYSIS_FILE);

        printGeneralStat();

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

    private static void printGeneralStat() throws FileNotFoundException {
        List<Set<String>> peptide = new ArrayList<>();
        for (int i = 0; i < Bioinf.peptide.size(); i++) {
            peptide.add(Bioinf.peptide.get(i).stream().filter(s -> !s.getSequences().isEmpty())
                    .flatMap(s -> s.getSequences().stream()).map(Sequence::getSequence).collect(Collectors.toSet()));
        }

        List<Set<String>> hcdPeptide = new ArrayList<>();
        for (int i = 0; i < Bioinf.peptide.size(); i++) {
            hcdPeptide.add(Bioinf.peptide.get(i).stream().filter(s -> !s.getSequences().isEmpty() && s.getActivation().equals(HCD))
                    .flatMap(s -> s.getSequences().stream()).map(Sequence::getSequence).collect(Collectors.toSet()));
        }

        List<Set<String>> cidPeptide = new ArrayList<>();
        for (int i = 0; i < Bioinf.peptide.size(); i++) {
            cidPeptide.add(Bioinf.peptide.get(i).stream().filter(s -> !s.getSequences().isEmpty() && s.getActivation().equals(CID))
                    .flatMap(s -> s.getSequences().stream()).map(Sequence::getSequence).collect(Collectors.toSet()));
        }

        List<Bag<Pair<Integer, String>>> heavyChainCoverages = new ArrayList<>();
        List<Bag<Pair<Integer, String>>> lightChainCoverages = new ArrayList<>();
        for (int i = 0; i < peptide.size(); i++) {
            Set<String> union = new HashSet<>();
            for (int j = 0; j < peptide.size(); j++) {
                if (i == j) {
                    continue;
                }
                union = Sets.union(union, peptide.get(j));
            }
            Set<String> diff = Sets.difference(hcdPeptide.get(i), union);
            heavyChainCoverages.add(getAppearancesFromSet(diff, heavyChain));
            lightChainCoverages.add(getAppearancesFromSet(diff, lightChain));
            diff = Sets.difference(cidPeptide.get(i), union);
            heavyChainCoverages.add(getAppearancesFromSet(diff, heavyChain));
            lightChainCoverages.add(getAppearancesFromSet(diff, lightChain));
        }

        try (PrintWriter writer = new PrintWriter(RESOURCES_FOLDER + GENERAL_STAT_FILE)) {
            for (int i = 0; i < lightChainCoverages.size(); i += 2) {
                System.out.print((i / 2 + 1) + " file, HCD: ");
                System.out.println(highlightCoverage(getCoverage(lightChainCoverages.get(i), lightChain.length()), lightChain));
                System.out.print((i / 2 + 1) + " file, CID: ");
                System.out.println(highlightCoverage(getCoverage(lightChainCoverages.get(i + 1), lightChain.length()), lightChain));

                writer.println((i / 2 + 1) + " file, HCD:");
                writer.println(lightChain);
                printFoundMatches(writer, lightChainCoverages.get(i));
                writer.println();
                writer.println((i / 2 + 1) + " file, CID:");
                writer.println(lightChain);
                printFoundMatches(writer, lightChainCoverages.get(i + 1));
                writer.println();
            }
            System.out.println();
            for (int i = 0; i < heavyChainCoverages.size(); i += 2) {
                System.out.print((i / 2 + 1) + " file, HCD: ");
                System.out.println(highlightCoverage(getCoverage(heavyChainCoverages.get(i), heavyChain.length()), heavyChain));
                System.out.print((i / 2 + 1) + " file, CID: ");
                System.out.println(highlightCoverage(getCoverage(heavyChainCoverages.get(i + 1), heavyChain.length()), heavyChain));

                writer.println((i / 2 + 1) + " file, HCD:");
                writer.println(heavyChain);
                printFoundMatches(writer, heavyChainCoverages.get(i));
                writer.println();
                writer.println((i / 2 + 1) + " file, CID:");
                writer.println(heavyChain);
                printFoundMatches(writer, heavyChainCoverages.get(i + 1));
                writer.println();
            }
            System.out.println();


            List<Set<String>> correctHcd20 = new ArrayList<>();
            for (int i = 0; i < Bioinf.peptide.size(); i++) {
                correctHcd20.add(Bioinf.peptide.get(i).stream().filter(s -> !s.getSequences().isEmpty() && s.getActivation().equals(HCD))
                        .flatMap(s -> s.getSequences().stream().limit(20)).map(Sequence::getSequence)
                        .filter(s -> !isSubstring(lightChain, s).isEmpty() || !isSubstring(heavyChain, s).isEmpty()).collect(Collectors.toSet()));
            }
            List<Set<String>> correctHcd2 = new ArrayList<>();
            for (int i = 0; i < Bioinf.peptide.size(); i++) {
                correctHcd2.add(Bioinf.peptide.get(i).stream().filter(s -> !s.getSequences().isEmpty() && s.getActivation().equals(HCD))
                        .flatMap(s -> s.getSequences().stream().limit(2)).map(Sequence::getSequence)
                        .filter(s -> !isSubstring(lightChain, s).isEmpty() || !isSubstring(heavyChain, s).isEmpty()).collect(Collectors.toSet()));
            }
            writer.println();
            printPeptidesLengthStat(writer, correctHcd20, correctHcd2);
        }
    }

    private static void printPeptidesLengthStat(@NotNull PrintWriter out, @NotNull List<Set<String>> peptide20, @NotNull List<Set<String>> peptide2) {
        List<Set<String>> uniquePeptide20 = getUnique(peptide20);
        List<Set<String>> uniquePeptide2 = getUnique(peptide2);
/*
        for (Set<Scan> aPeptide : peptide) {
            peptide20.add(aPeptide.stream().filter(s -> s.getActivation().equals(HCD))
                    .flatMap(s -> s.getSequences().stream().limit(20)).map(Sequence::getSequence).collect(Collectors.toList()));
            peptide2.add(aPeptide.stream().filter(s -> s.getActivation().equals(HCD))
                    .flatMap(s -> s.getSequences().stream().limit(2)).map(Sequence::getSequence).collect(Collectors.toList()));
        }
*/
        List<String> longest2 = new ArrayList<>();
        List<String> longest20 = new ArrayList<>();
        for (int i = 0; i < peptide.size(); i++) {
            int j = i;
            out.println((i + 1) + " file, 20 interpretations, longest:");
            uniquePeptide20.get(i).stream().sorted(Comparator.comparingInt(String::length).reversed()).limit(3).peek(out::println).filter(s -> j == 0).forEach(longest20::add);
            out.println((i + 1) + " file, 20 interpretations, shortest:");
            uniquePeptide20.get(i).stream().sorted(Comparator.comparingInt(String::length)).limit(3).forEach(out::println);
            out.print((i + 1) + " file, 20 interpretations, average: ");
            out.println(uniquePeptide20.get(i).stream().mapToInt(String::length).average().orElse(0));
            out.println();

            out.println((i + 1) + " file, 2 interpretations, longest:");
            uniquePeptide2.get(i).stream().sorted(Comparator.comparingInt(String::length).reversed()).limit(3).peek(out::println).filter(s -> j == 0).forEach(longest2::add);
            out.println((i + 1) + " file, 2 interpretations, shortest:");
            uniquePeptide2.get(i).stream().sorted(Comparator.comparingInt(String::length)).limit(3).forEach(out::println);
            out.print((i + 1) + " file, 2 interpretations, average: ");
            out.println(uniquePeptide2.get(i).stream().mapToInt(String::length).average().orElse(0));
            out.println();
        }

        for(String longest : longest2) {
            out.println("2 interpretations, 2 file coverage:");
            coverWithPeptideFormSet(out, longest, uniquePeptide2.get(1));
            out.println("2 interpretations, 3 file coverage:");
            coverWithPeptideFormSet(out, longest, uniquePeptide2.get(2));
        }

        for(String longest : longest20) {
            out.println("20 interpretations, 2 file coverage:");
            coverWithPeptideFormSet(out, longest, uniquePeptide20.get(1));
            out.println("20 interpretations, 3 file coverage:");
            coverWithPeptideFormSet(out, longest, uniquePeptide20.get(2));
        }
    }

    private static void coverWithPeptideFormSet(PrintWriter out, String longest, Set<String> strings) {
        List<Integer> lightPos = isSubstring(lightChain, longest);
        List<Integer> heavyPos = isSubstring(heavyChain, longest);

        Bag<Pair<Integer, String>> lightCover = new TreeBag<>(Comparator.comparing(Pair::getKey));
        Bag<Pair<Integer, String>> heavyCover = new TreeBag<>(Comparator.comparing(Pair::getKey));
        for(String sequence : strings) {
            List<Integer> positions = isSubstring(lightChain, sequence);
            if (isIntersect(lightPos, longest, positions, sequence)) {
                for (int position : positions) {
                    lightCover.add(new Pair<>(position, sequence));
                }
            }

            positions = isSubstring(heavyChain, sequence);
            if (isIntersect(heavyPos, longest, positions, sequence)) {
                for (int position : positions) {
                    heavyCover.add(new Pair<>(position, sequence));
                }
            }
        }

        if (lightPos.size() != 0) {
            out.println(lightChain);
            printWhitespaces(out, lightPos.get(0));
            out.println(longest);
            printFoundMatches(out, lightCover);
        }

        if (heavyPos.size() != 0) {
            out.println(heavyChain);
            printWhitespaces(out, heavyPos.get(0));
            out.println(longest);
            printFoundMatches(out, heavyCover);
        }

        out.println();
    }

    private static boolean isIntersect(List<Integer> pos1, String seq1, List<Integer> pos2, String seq2) {
        boolean res = false;

        for(int i : pos1) {
            for(int j : pos2) {
                res |= i <= j && i + seq1.length() > j;
                res |= i < j + seq2.length() && i + seq1.length() >= j + seq2.length();
                res |= i >= j && i + seq1.length() <= j + seq2.length();
                }
        }

        return res;
    }


    private static List<Set<String>> getUnique(List<Set<String>> notUnique) {
        List<Set<String>> res = new ArrayList<>();
        for (int i = 0; i < notUnique.size(); i++) {
            Set<String> union = new HashSet<>();
            for (int j = 0; j < notUnique.size(); j++) {
                if (i == j) {
                    continue;
                }
                union = Sets.union(union, notUnique.get(j));
            }
            res.add(Sets.difference(notUnique.get(i), union));
        }

        return res;
    }

    private static Bag<Pair<Integer, String>> getAppearancesFromSet(Set<String> set, String chain) {
        Bag<Pair<Integer, String>> appearances = new TreeBag<>(Comparator.comparing(Pair::getKey));

        for (String sequence : set) {
            List<Integer> positions = isSubstring(chain, sequence);
            if (!positions.isEmpty()) {
                for (int position : positions) {
                    appearances.add(new Pair<>(position, sequence));
                }
            }
        }

        return appearances;/*
        BitSet coverage = getCoverage(appearances, chain.length());
        return highlightCoverage(coverage, chain);*/
    }

    private static void savePeptide() {
        peptide.add(new HashSet<>(scans.values()));
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
            if (scan.getActivation().equals(CID)) {
                cidScans++;
            } else {
                hcdScans++;
            }
            for (int i = 0; i < Math.min(numberOfInterpretations + 1, scan.getSequences().size()); i++) {
                Sequence sequence = scan.getSequences().get(i);
                if (!isSubstring(lightChain, sequence.getSequence()).isEmpty()) {
                    if (scan.getActivation().equals(HCD)) {
                        hcdCorrectInterpretations++;
                    } else {
                        cidCorrectInterpretations++;
                    }
                } else if (!isSubstring(heavyChain, sequence.getSequence()).isEmpty()) {
                    if (scan.getActivation().equals(HCD)) {
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

    @NotNull
    private static Bag<Pair<Integer, String>> getAppearances(@NotNull String chain, @NotNull Activation activation, int numberOfTries) {
        Bag<Pair<Integer, String>> appearances = new TreeBag<>(Comparator.comparing(Pair::getKey));

        for (Scan scan : scans.values()) {
            for (int i = 0; i < Math.min(numberOfTries + 1, scan.getSequences().size()); i++) {
                Sequence sequence = scan.getSequences().get(i);
                List<Integer> positions = isSubstring(chain, sequence.getSequence());
                if (!positions.isEmpty() && scan.getActivation().equals(activation)) {
                    for (int position : positions) {
                        appearances.add(new Pair<>(position, sequence.getSequence()));
                    }
                }
            }
        }

        return appearances;
    }

    @NotNull
    private static Coverage printMatches(PrintWriter out, int numberOfCidTries, int numberOfHcdTries) {
        Bag<Pair<Integer, String>> hcdLightChain = getAppearances(lightChain, HCD, numberOfHcdTries);
        Bag<Pair<Integer, String>> hcdHeavyChain = getAppearances(heavyChain, HCD, numberOfHcdTries);
        Bag<Pair<Integer, String>> cidLightChain = getAppearances(lightChain, CID, numberOfCidTries);
        Bag<Pair<Integer, String>> cidHeavyChain = getAppearances(heavyChain, CID, numberOfCidTries);

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

    @NotNull
    private static String highlightCoverage(BitSet coverage, String chain) {
        StringBuilder highlighted = new StringBuilder();

        for (int i = 0; i < chain.length(); i++) {
            if (coverage.get(i)) {
                highlighted.append(ConsoleColors.GREEN).append(chain.charAt(i)).append(ConsoleColors.RESET);
            } else {
                highlighted.append(chain.charAt(i));
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
        if (pair != null && Math.abs(pair.getPrecursorMass() - mass) <= EPSILON) {
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

    @NotNull
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

