package nl.systemsgenetics.gwassummarystatistics;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import com.opencsv.CSVWriter;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.Allele;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class MatrixBasedGwasSummaryStatistics implements GwasSummaryStatistics{

    private static final Logger LOGGER = Logger.getLogger(WritableGwasSummaryStatistics.class);
    private String gwasId;
    private final LinkedHashMap<String, Integer> variantIds;
    private final DoubleMatrix1D betaCoefficients;
    private final DoubleMatrix1D pValues;

    public MatrixBasedGwasSummaryStatistics(String gwasId, LinkedHashMap<String, Integer> variantIds,
                                            DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues) {
        this.gwasId = gwasId;
        this.variantIds = variantIds;
        this.betaCoefficients = betaCoefficients;
        this.pValues = pValues;
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    @Override
    public float[] getEffectSizeEstimates(GeneticVariant variant) {
        return new float[]{(float) betaCoefficients.getQuick(variantIds.get(variant.getPrimaryVariantId()))};
    }

    @Override
    public float[] getTransformedPValues(GeneticVariant variant) {
        return new float[]{(float) -Math.log10(pValues.getQuick(variantIds.get(variant.getPrimaryVariantId())))};
    }

    @Override
    public Iterator<EffectAllele> effectAlleles() {
        throw new UnsupportedOperationException("Not currently supported");
    }

    @Override
    public THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> riskEntries(
            RandomAccessGenotypeData genotypeData, double[] pValThres,
            String[] genomicRangesToExclude, boolean unweighted) {
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
        THashMap<String, THashMap<String, THashMap<String, ArrayList<RiskEntry>>>> risks = new THashMap<>();

        THashSet<String> chromosomesExcluded = new THashSet<>();
        int snpsExcluded = 0;

        THashMap<String, THashMap<String, ArrayList<RiskEntry>>> filehash = new THashMap<String, THashMap<String, ArrayList<RiskEntry>>>();
        String name = getGwasId();

        for (double p : pValThres) {
            String name2 = "_P" + p;
            if (!filehash.containsKey(name2)) {
                filehash.put(name2, new THashMap<>());
            }
        }

        for (Map.Entry<String, Integer> variantAssociation : variantIds.entrySet()) {

//                    System.out.println(s);
            if (genotypeData.getVariantIdMap().containsKey(variantAssociation.getKey())) {
                GeneticVariant snpObject = genotypeData.getVariantIdMap().get(variantAssociation.getKey());
//                        System.out.print(snpObject.getSequenceName() + "\t" + snpObject.getStartPos() + "\n");
                double currentP = pValues.getQuick(variantAssociation.getValue());
                boolean addEntry = true;

                double or = betaCoefficients.getQuick(variantAssociation.getValue());

                if (unweighted) {
                    if (or < 0) {
                        or = -1;
                    } else {
                        or = 1;
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
                                filehash.get(name2).put(snpObject.getSequenceName(), new ArrayList<>());
                            }
                            filehash.get(name2).get(snpObject.getSequenceName()).add(new RiskEntry(
                                    variantAssociation.getKey(), snpObject.getSequenceName(), snpObject.getStartPos(),
                                    snpObject.getAlternativeAlleles()
                                            .get(snpObject.getAlternativeAlleles().getAlleleCount() - 1).getAlleleAsSnp(),
                                    or, currentP));
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

    public void save(String pathPrefix) throws IOException {
        // Format the output path.
        String outputPath = String.format("%s_%s.tsv",
                pathPrefix, gwasId.replace(" ", "_"));

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Write a header (is skipped when reading again)
        writer.writeNext(new String[]{"variant", "AssessedAllele", "ES", "pvalue"});

        // Write risk entries
        for (Map.Entry<String, Integer> variantAssociation : variantIds.entrySet()) {
            writer.writeNext(new String[]{
                    variantAssociation.getKey(),
                    "LastAlternativeAllele",
                    String.valueOf(betaCoefficients.getQuick(variantAssociation.getValue())),
                    String.valueOf(pValues.getQuick(variantAssociation.getValue()))});
        }

        writer.close();
    }

    public int size() {
        return variantIds.size();
    }
}
