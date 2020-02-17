package nl.systemsgenetics.gwassummarystatistics;

import com.google.common.primitives.Floats;
import com.opencsv.CSVWriter;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.polygenicscorecalculator.Main;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class ReadOnlyLegacyGwasSummaryStatistics implements GwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(ReadOnlyLegacyGwasSummaryStatistics.class);
    private static final Pattern TAB_PATTERN = Pattern.compile("\\t");
    private String gwasId;
    private String riskFilePath;

    public ReadOnlyLegacyGwasSummaryStatistics(String gwasId, String riskFilePath) {
        this.gwasId = gwasId;
        this.riskFilePath = riskFilePath;
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    public float[] getEffectSizeEstimates(GeneticVariant variant) {
        throw new UnsupportedOperationException("Not supported");
    }

    public float[] getTransformedPValues(GeneticVariant variant) {
        throw new UnsupportedOperationException("Not supported");
    }

    public Iterator<EffectAllele> effectAlleles() {
        throw new UnsupportedOperationException("Not currently supported");
    }

    @Override
    public THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            RandomAccessGenotypeData genotypeData,
            double[] pValThres,
            String[] genomicRangesToExclude,
            boolean unweighted) {
        THashMap<String, ArrayList<Pair<Integer, Integer>>> exclussionRanges = new THashMap<>();

        if (genomicRangesToExclude != null) {
            System.out.println("Trying to exclude genomic ranges.");
            int ranges = 0;
            for (String s : genomicRangesToExclude) {
                String[] parts = s.split(":");
                String key = parts[0];
                parts = parts[1].split("-");

                if (!exclussionRanges.contains(key)) {
                    exclussionRanges.put(key, new ArrayList<>());
                }
                exclussionRanges.get(key).add(new Pair(Integer.parseInt(parts[0]), Integer.parseInt(parts[1])));
                ranges++;
            }
            if (LOGGER.isDebugEnabled()) {
                System.out.println("Number of ranges excluded: " + ranges + " on: " + exclussionRanges.size() + " chromosomes");
            }
        }
        // A risk entry (value?) per variant? per sequence? per pval threshold? per file?
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = new THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>>();

        File riskFile = new File(this.riskFilePath);
        if (!riskFile.exists()) {
            System.out.println("Warning: input risk folder does not exists:\n");
            System.out.println(riskFile);
            System.exit(-1);
        }

        THashSet<String> chromosomesExcluded = new THashSet<>();
        int snpsExcluded = 0;
        try {
            TextFile readFiles = new TextFile(riskFile.getAbsolutePath(), TextFile.R);

            System.out.println(riskFile.getName());


            THashMap<String, THashMap<String, ArrayList<RiskEntry>>> filehash = new THashMap<String, THashMap<String, ArrayList<RiskEntry>>>();
            String name = riskFile.getName();

            for (double p : pValThres) {
                String name2 = "_P" + p;
                if (!filehash.containsKey(name2)) {
                    filehash.put(name2, new THashMap<String, ArrayList<RiskEntry>>());
                }
            }

//				synchronized (risks) {
//					risks.put(name, filehash);
//					for (double p : pValueThreshold) {
//
//						String name2 = "_P" + p;
//
//						if (!risks.containsKey(name)) {
//							risks.put(name, new THashMap<String, THashMap<String, ArrayList<RiskEntry>>>());
//						}
//
//						if (!risks.get(name).containsKey(name2)) {
//							risks.get(name).put(name2, new THashMap<String, ArrayList<RiskEntry>>());
//						}
//					}
//				}

            String s = readFiles.readLine();
            while ((s = readFiles.readLine()) != null) {
                String[] parts = TAB_PATTERN.split(s);
//                    System.out.println(s);
                if (genotypeData.getVariantIdMap().containsKey(parts[0])) {
                    GeneticVariant snpObject = genotypeData.getVariantIdMap().get(parts[0]);
//                        System.out.print(snpObject.getSequenceName() + "\t" + snpObject.getStartPos() + "\n");
                    double currentP = Double.parseDouble(parts[3]);
                    boolean addEntry = true;

                    if (unweighted) {
                        if (parts[2].startsWith("-")) {
                            parts[2] = "-1";
                        } else {
                            parts[2] = "1";
                        }
                    }

                    if (exclussionRanges.contains(snpObject.getSequenceName())) {
                        chromosomesExcluded.add(snpObject.getSequenceName());
                        for (Pair<Integer, Integer> p : exclussionRanges.get(snpObject.getSequenceName())) {
                            if (p.getLeft() <= snpObject.getStartPos() && p.getRight() >= snpObject.getStartPos()) {
                                addEntry = false;
                                snpsExcluded++;
                            }
                        }
                    }

                    if (addEntry) {
                        for (double p : pValThres) {
                            if (currentP < p) {
                                String name2 = "_P" + p;

                                if (!filehash.get(name2).containsKey(snpObject.getSequenceName())) {
                                    filehash.get(name2).put(snpObject.getSequenceName(), new ArrayList<RiskEntry>());
                                }
                                filehash.get(name2).get(snpObject.getSequenceName()).add(new RiskEntry(parts[0],
                                        snpObject.getSequenceName(), snpObject.getStartPos(),
                                        parts[1].charAt(0), Double.parseDouble(parts[2]), currentP));
                            }
                        }
                    }
                }
            }
            synchronized (risks) {
                risks.put(name, filehash);
            }

            if (LOGGER.isDebugEnabled()) {
                System.out.println("Chromosomes where regions are excluded: " + chromosomesExcluded);
                System.out.println("Number of SNPs excluded: " + snpsExcluded);
            }
            readFiles.close();
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }

        for (Map.Entry<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> e : risks.entrySet()) {

            for (Map.Entry<String, THashMap<String, ArrayList<RiskEntry>>> e2 : e.getValue().entrySet()) {
                int entries = 0;
                for (Map.Entry<String, ArrayList<RiskEntry>> e3 : e2.getValue().entrySet()) {
                    Collections.sort(e3.getValue());
                    entries += e3.getValue().size();
                }
//                System.out.println(e.getKey()+e2.getKey()+" has: "+entries+" entries");
            }

        }

        return risks;
    }
}
