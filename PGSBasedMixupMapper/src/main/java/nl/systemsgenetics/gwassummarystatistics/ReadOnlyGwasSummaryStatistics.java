package nl.systemsgenetics.gwassummarystatistics;

import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;

import java.util.*;

public class ReadOnlyGwasSummaryStatistics implements GwasSummaryStatistics {

    private static final Logger LOGGER = Logger.getLogger(ReadOnlyGwasSummaryStatistics.class);
    private MultiStudyGwasSummaryStatistics originalSummaryStatistics;
    private int studyIndex;
    private String studyName;

    public ReadOnlyGwasSummaryStatistics(MultiStudyGwasSummaryStatistics originalSummaryStatistics, String studyName) {

        List<String> sampleNames = Arrays.asList(originalSummaryStatistics.getSampleNames());
        if (!sampleNames.contains(studyName)) {
            throw new GwasSummaryStatisticsException(
                    String.format("Study name %s does not exist in the given summary statistics", studyName));
        }

        this.originalSummaryStatistics = originalSummaryStatistics;
        this.studyName = studyName;
        this.studyIndex = sampleNames.indexOf(studyName);
    }

    @Override
    public String getGwasId() {
        return studyName;
    }

    @Override
    public float[] getEffectSizeEstimates(GeneticVariant variant) {
        return originalSummaryStatistics.getEffectSizeEstimates(variant)[studyIndex];
    }

    @Override
    public float[] getTransformedPValues(GeneticVariant variant) {
        return originalSummaryStatistics.getTransformedPValues(variant)[studyIndex];
    }

    @Override
    public Iterator<EffectAllele> effectAlleles() {
        return originalSummaryStatistics.effectAlleles(this);
    }

    @Override
    public THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            RandomAccessGenotypeData genotypeData,
            double[] pValueThreshold, String[] genomicRangesToExclude, boolean unweighted) {
        THashMap<String, ArrayList<Pair<Integer, Integer>>> exclussionRanges = new THashMap<>();

        if (genomicRangesToExclude != null) {
            System.out.println("Trying to exclude genomic ranges.");
            int ranges = 0;
            for (String s : genomicRangesToExclude) {
                String[] parts = s.split(":");
                String key = parts[0];
                parts = parts[1].split("-");

                if (!exclussionRanges.contains(key)) {
                    exclussionRanges.put(key, new ArrayList<Pair<Integer, Integer>>());
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

        THashSet<String> chromosomesExcluded = new THashSet<>();
        int snpsExcluded = 0;
        try {
            String name = this.getGwasId();

            THashMap<String, THashMap<String, ArrayList<RiskEntry>>> filehash = new THashMap<>();

            for (double p : pValueThreshold) {
                String name2 = "_P" + p;
                if (!filehash.containsKey(name2)) {
                    filehash.put(name2, new THashMap<>());
                }
            }

            int numberOfVariantsWithoutSameAlleles = 0;
            int numberOfVariantsWithoutSameAllelesComplement = 0;

            for (GeneticVariant variant : originalSummaryStatistics) {
                // Check if there is a variant in the genotype data corresponding to this variant
                GeneticVariant snpVariantByPos = genotypeData.getSnpVariantByPos(
                        variant.getSequenceName(),
                        variant.getStartPos());

                if (snpVariantByPos == null ||
                        !snpVariantByPos.getPrimaryVariantId().equals(variant.getPrimaryVariantId())) {
                    continue;
                }

                if (!snpVariantByPos.getVariantAlleles().sameAlleles(variant.getVariantAlleles())) {
                    numberOfVariantsWithoutSameAlleles++;
                    Alleles complement = snpVariantByPos.getVariantAlleles().getComplement();

                    if (!complement.sameAlleles(variant.getVariantAlleles())) {
                        numberOfVariantsWithoutSameAllelesComplement++;
                        continue;
                    }
                }

//                    System.out.println(s);
//                        System.out.print(snpObject.getSequenceName() + "\t" + snpObject.getStartPos() + "\n");
                double currentP = Math.pow(10, -this.getTransformedPValues(variant)[0]);
                boolean addEntry = true;
                float partsTwo = this.getEffectSizeEstimates(variant)[0];

                if (unweighted) {
                    if (partsTwo < 0) {
                        partsTwo = -1;
                    } else {
                        partsTwo = 1;
                    }
                }

                if (exclussionRanges.contains(variant.getSequenceName())) {
                    chromosomesExcluded.add(variant.getSequenceName());
                    for (Pair<Integer, Integer> p : exclussionRanges.get(variant.getSequenceName())) {
                        if (p.getLeft() <= variant.getStartPos() && p.getRight() >= variant.getStartPos()) {
                            addEntry = false;
                            snpsExcluded++;
                        }
                    }
                }

                if (addEntry) {
                    for (double p : pValueThreshold) {
                        if (currentP < p) {
                            String name2 = "_P" + p;

                            if (!filehash.get(name2).containsKey(variant.getSequenceName())) {
                                filehash.get(name2).put(variant.getSequenceName(), new ArrayList<>());
                            }
                            filehash.get(name2).get(variant.getSequenceName()).add(new RiskEntry(variant.getPrimaryVariantId(),
                                    variant.getSequenceName(), variant.getStartPos(),
                                    variant.getAlternativeAlleles().getAllelesAsChars()[0], partsTwo, currentP));
                        }
                    }
                }
            }
            synchronized (risks) {
                risks.put(name, filehash);
            }

            if (numberOfVariantsWithoutSameAlleles > 0) {
                String message = String.format(
                        "%d variants have different alleles than the alleles of " +
                                "the matched variant from the input genotype data.%n" +
                                "The complement of these alleles is also different for %d variants, " +
                                "which will be removed.",
                        numberOfVariantsWithoutSameAlleles, numberOfVariantsWithoutSameAllelesComplement);
                System.err.println(message);
                LOGGER.warn(message);
            }

            if (LOGGER.isDebugEnabled()) {
                System.out.println("Chromosomes where regions are excluded: " + chromosomesExcluded);
                System.out.println("Number of SNPs excluded: " + snpsExcluded);
            }
        } catch (GwasSummaryStatisticsException ex) {
            LOGGER.error(ex);
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
