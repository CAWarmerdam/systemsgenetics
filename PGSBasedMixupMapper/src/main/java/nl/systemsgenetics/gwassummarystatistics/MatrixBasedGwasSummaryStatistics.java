package nl.systemsgenetics.gwassummarystatistics;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import com.opencsv.CSVWriter;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import nl.systemsgenetics.polygenicscorecalculator.RiskEntry;
import org.apache.log4j.Logger;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import umcg.genetica.containers.Pair;

import javax.annotation.Nonnull;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class MatrixBasedGwasSummaryStatistics implements GwasSummaryStatistics{

    private static final Logger LOGGER = Logger.getLogger(MatrixBasedGwasSummaryStatistics.class);
    private String gwasId;
    private final LinkedHashMap<String, Integer> variantIds;
    private DoubleMatrix1D betaCoefficients;
    private DoubleMatrix1D pValues;

    public MatrixBasedGwasSummaryStatistics(String gwasId, LinkedHashMap<String, Integer> variantIds,
                                            DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues) {
        // Set the gwas identifier
        this.gwasId = gwasId;

        // Make sure that the sizes of the three collections are equal.
        if (variantIds.size() != betaCoefficients.size() || variantIds.size() != pValues.size()) {
            throw new GwasSummaryStatisticsException(
                    String.format("Trying to add summary statistics with unequal number of variants, " +
                                    "betas or p-values (%d vs %d, %d)",
                            variantIds.size(), betaCoefficients.size(), pValues.size()));
        }

        this.variantIds = variantIds;
        this.betaCoefficients = betaCoefficients;
        this.pValues = pValues;
    }

    public MatrixBasedGwasSummaryStatistics(RandomAccessGenotypeData genotypeData, String gwasId) {
        this.gwasId = gwasId;
        this.variantIds = new LinkedHashMap<>();
        this.betaCoefficients = new DenseDoubleMatrix1D(0);
        this.pValues = new DenseDoubleMatrix1D(0);
    }

    public void add(LinkedHashMap<String, Integer> variantIds, DoubleMatrix1D betaCoefficients, DoubleMatrix1D pValues) {
        if (variantIds.size() != betaCoefficients.size() || variantIds.size() != pValues.size()) {
            throw new GwasSummaryStatisticsException(
                    String.format("Trying to add summary statistics with unequal number of variants, " +
                            "betas or p-values (%d vs %d, %d)",
                            variantIds.size(), betaCoefficients.size(), pValues.size()));
        }

        this.betaCoefficients = DoubleFactory1D.dense.append(this.betaCoefficients, betaCoefficients);
        this.pValues = DoubleFactory1D.dense.append(this.pValues, pValues);

        int originalSize = this.variantIds.size();

        for (String variantId : variantIds.keySet()) {
            this.variantIds.put(variantId, variantIds.get(variantId) + originalSize);
        }
    }

    @Override
    public String getGwasId() {
        return gwasId;
    }

    @Override
    @Nonnull
    public Iterator<EffectAllele> iterator() {
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

    public void save(Path pathPrefix, Double pValueThreshold) throws IOException {
        // Format the output path.
        Path outputPath = Paths.get(String.format("%s.txt",
                pathPrefix));

        LOGGER.info(String.format("Writing summary statistics to '%s'", outputPath));

        if (!outputPath.getParent().toFile().isDirectory() && !outputPath.getParent().toFile().mkdir()) {
            throw new IOException(String.format("Could not create directory '%s'", outputPath.getParent()));
        }

        // Create a CSV writer without quotation characters and with tabs as separators to mimic
        // the legacy gwas summary statistics files.
        CSVWriter writer = new CSVWriter(new FileWriter(outputPath.toFile()),
                '\t', CSVWriter.NO_QUOTE_CHARACTER,
                CSVWriter.DEFAULT_ESCAPE_CHARACTER, CSVWriter.DEFAULT_LINE_END);

        // Write a header (is skipped when reading again)
        writer.writeNext(new String[]{"variant", "AssessedAllele", "ES", "pvalue"});

        // Write risk entries
        for (Map.Entry<String, Integer> variantAssociation : variantIds.entrySet()) {

            // Check if the pvalue for this association is below the given threshold.
            if (pValues.getQuick(variantAssociation.getValue()) < pValueThreshold) {
                // Only write the association whenever the p-value is below the threshold
                writer.writeNext(new String[]{
                        variantAssociation.getKey(),
                        "LastAlternativeAllele",
                        String.valueOf(betaCoefficients.getQuick(variantAssociation.getValue())),
                        String.valueOf(pValues.getQuick(variantAssociation.getValue()))});
            }
        }

        writer.close();
    }

    public int size() {
        return variantIds.size();
    }
}
